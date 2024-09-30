from flask import Flask, render_template, request, jsonify
from multiprocessing import Process, Manager, Queue
from lverage import Lverage, LverageRecord
from src.MotifDB import MotifDBFactory
import os
import src.exceptions as lvexceptions

import pandas as pd
from collections import OrderedDict

app = Flask(__name__)

manager = Manager()
lverage_dict = manager.dict()
queue = Queue()

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/get_motif_db_names', methods=['GET'])
def get_motif_db_names():
    db_names = MotifDBFactory.get_motif_db_names()
    return jsonify(db_names)

@app.route('/run', methods=['POST'])
def run():
    data = request.json
    email = data.get('email')
    orthologs = data.get('orthologs')
    motif_db = data.get('motif_db')
    project_title = data.get('project_title')
    fasta_directory = data.get('fasta_directory')
    output_directory = data.get('output_directory')

    error_messages = []

    if project_title == '':
        error_messages.append({'field': 'project-title', 'message': 'Project title cannot be empty'})

    if not os.path.isdir(fasta_directory):
        error_messages.append({'field': 'fasta-directory', 'message': 'FASTA directory does not exist'})
    else:
        valid_extensions = ['.fa', '.fasta', '.fna']
        fasta_files = [os.path.join(fasta_directory, f) for f in os.listdir(fasta_directory) 
                            if os.path.isfile(os.path.join(fasta_directory, f)) 
                            and os.path.splitext(f)[1] in valid_extensions]
        fasta_files.sort(key=lambda x: os.path.basename(x))

        if len(fasta_files) == 0:
            error_messages.append({'field': 'fasta-directory', 'message': 'No FASTA files found in directory.'})

    if not os.path.isdir(output_directory):
        error_messages.append({'field': 'output-directory', 'message': 'Output directory does not exist'})
    
    exception_list = Lverage.check_valid(motif_database=motif_db, ortholog_name_list=orthologs, email=email)
    for exception in exception_list:
        if isinstance(exception, lvexceptions.OrthologSpeciesError):
            error_messages.append({'field': 'ortholog-species', 'message': str(exception)})
        elif isinstance(exception, lvexceptions.EmailError):
            error_messages.append({'field': 'email', 'message': str(exception)})
        else:
            error_messages.append({'field': 'run-button', 'message': str(exception)})

    output_file = os.path.join(output_directory, project_title.strip() + '_output.tsv')

    if os.path.exists(output_file):
        error_messages.append({'field': 'project-title', 'message': f'Project file already exists: {output_file}'})

    if len(error_messages) > 0:
        return jsonify({'status': 'error', 'errors': error_messages})

    lverage = Lverage(motif_database=motif_db, ortholog_name_list=orthologs, email=email, verbose=True)

    lverage_dict['lverage'] = lverage
    lverage_dict['output_file'] = output_file
    lverage_dict['cancel'] = False

    f = open(output_file, 'w')
    headers = ["Gene ID"] + LverageRecord.header_list + MotifDBFactory.get_db_record(motif_db).header_list
    f.write('\t'.join(headers) + '\n')
    f.close()

    result = {'status': 'progress', 'fasta_files': fasta_files, 'project_title': project_title}
    return jsonify(result)

def process_file_task(fasta_file, lverage_dict, queue):
    lverage = lverage_dict['lverage']
    gene_id = os.path.basename(fasta_file).split('.')[0]

    try:
        motifs = lverage.call(gene_file=fasta_file)
        diagnostics = lverage.diagnostic_list

        for motif, diagnostic in zip(motifs, diagnostics):
            if lverage_dict.get('cancel'):
                queue.put({'status': 'cancelled'})
                return

            motif_values = motif.get_values()
            diagnostic_values = diagnostic.get_values()
            motif_values = [str(x) for x in motif_values]
            diagnostic_values = [str(x) for x in diagnostic_values]

            values = [gene_id] + diagnostic_values + motif_values
            queue.put({'status': 'success', 'data': '\t'.join(values) + '\n'})
        
    except Exception as e:
        queue.put({'status': 'error', 'error': str(e)})

@app.route('/process_file', methods=['POST'])
def process_file():
    data = request.json
    fasta_file = data.get('fasta_file')
    queue = Queue()

    lverage_dict['cancel'] = False

    p = Process(target=process_file_task, args=(fasta_file, lverage_dict, queue))
    p.start()
    
    while p.is_alive():
        p.join(timeout=0.1)  # Adjust the timeout as necessary
        if lverage_dict['cancel']:
            p.terminate()
            p.join()  # Ensure the process has terminated
    output_file = lverage_dict['output_file']

    while not queue.empty():
        result = queue.get()
        if result['status'] == 'error':
            return jsonify({'status': 'error', 'error': result['error']})
        elif result['status'] == 'cancelled':
            return jsonify({'status': 'cancelled'})
        elif result['status'] == 'success':
            with open(output_file, 'a') as f:
                f.write(result['data'])
                f.flush()

    return jsonify({'status': 'success'})

@app.route('/cancel_run', methods=['POST'])
def cancel_run():
    lverage_dict['cancel'] = True
    return jsonify({'status': 'success'})

@app.route('/close_session', methods=['POST'])
def close_session():
    if 'output_file' in lverage_dict and lverage_dict['output_file']:
        with open(lverage_dict['output_file'], 'a') as f:
            f.close()

    lverage_dict['output_file'] = None
    lverage_dict['lverage'] = None
    lverage_dict['cancel'] = False

    return jsonify({'status': 'success'})

@app.route('/progress')
def progress():
    return render_template('progress.html')

@app.route('/validate_output', methods=['POST'])
def validate_output():
    data = request.json
    output_content = data.get('file_content')

    lines = [line for line in output_content.split('\n') if line != '']

    # Check if there is at least one line (not including header)


    if len(lines) <= 1:
        return jsonify({'status': 'error', 'message': 'Output file is empty'})
 
    # all good
    return jsonify({'status': 'success'})

@app.route('/prepare_output', methods=['POST'])
def prepare_output():
    data = request.json
    output_content = data.get('file_content')

    try:
        lines = output_content.split('\n')
        headers = lines[0].split('\t')

        lverage_dict['headers'] = headers

        data = [line.split('\t') for line in lines[1:] if line != '']

        df = pd.DataFrame(data, columns=headers)

        motif_dict = OrderedDict()

        grouped = df.groupby(['Gene ID', 'Gene DBD Name', 'Ortholog BLAST Description'])
        for (gene_id, dbd, description), group in grouped:
            if gene_id not in motif_dict:
                motif_dict[gene_id] = OrderedDict()
            if dbd not in motif_dict[gene_id]:
                motif_dict[gene_id][dbd] = OrderedDict()
            if description not in motif_dict[gene_id][dbd]:
                motif_dict[gene_id][dbd][description] = []
            motif_dict[gene_id][dbd][description].extend(group.to_dict('records'))

        # Sort the motif_dict by Gene ID, Gene DBD Name, and BLAST Description
        motif_dict = OrderedDict(sorted(motif_dict.items()))
        for gene_id in motif_dict:
            motif_dict[gene_id] = OrderedDict(sorted(motif_dict[gene_id].items()))
            for dbd in motif_dict[gene_id]:
                motif_dict[gene_id][dbd] = OrderedDict(sorted(motif_dict[gene_id][dbd].items()))


        return jsonify({'status': 'success', 'motif_dict': motif_dict})
    
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)})

@app.route('/get_motifs_info', methods=['POST'])
def get_motifs_info():
    data = request.json
    lines = data.get('motif_lines')


    motif_info = []
    motif_links = []
    for record in lines:

        motif_links.append(record['Motif Logo Link'])

        s = ""
        for key, value in record.items():
            s += f"{key}: {value}\n"
        motif_info.append(s)

    return jsonify({'status': 'success', 'motif_info': motif_info, 'motif_links': motif_links})


@app.route('/output')
def output():
    return render_template('output.html')

if __name__ == '__main__':
    app.run(debug=False)
