#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
The purpose of this file is to provide a list of exceptions that can be raised in the program.

Copyright (C) <RELEASE_YEAR_HERE> Bradham Lab

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

Correspondence: anthonygarza124@gmail.com OR ...cyndi's email here...
'''
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

##@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
from enum import Enum

##@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Exception Classes

class LverageErrorCode(Enum):
    '''
    This class contains the error codes that can be raised in the program
    '''
    MOTIF_DATABASE_NOT_FOUND = "Motif database not found"
    ORTHOLOG_SPECIES_NOT_FOUND = "Ortholog species not found in NCBI database"
    INVALID_EMAIL = "Invalid email address"
    INVALID_IDENTITY_THRESHOLD = "Invalid identity threshold. Must be between 0 and 1"
    CANT_CONVERT_TO_TAX_ID = "Can't convert ortholog to taxon ID"
    INVALID_BLAST_HIT_COUNT = "Invalid number of hits for BlastP. Must be greater than 0."
    INVALID_MOTIF_HIT_COUNT = "Invalid number of hits for motif database. Must be greater than 0."
    INVALID_ESCORE_THRESHOLD = "Invalid e-score threshold. Must be greater than 0."

class LverageError(Exception):
    '''
    This class is the parent class for all exceptions in the program
    '''

    def __init__(self, message, code: LverageErrorCode):
        super().__init__(message)
        self.code = code

class MotifDatabaseError(LverageError):
    '''
    Exception raised for errors in the motif database.
    '''

    def __init__(self, message=LverageErrorCode.MOTIF_DATABASE_NOT_FOUND.value):
        super().__init__(message, LverageErrorCode.MOTIF_DATABASE_NOT_FOUND)

class OrthologSpeciesError(LverageError):
    '''
    Exception raised for errors in ortholog species.
    '''

    def __init__(self, message=LverageErrorCode.ORTHOLOG_SPECIES_NOT_FOUND.value):
        super().__init__(message, LverageErrorCode.ORTHOLOG_SPECIES_NOT_FOUND)


class EmailError(LverageError):
    '''
    Exception raised for invalid email address.
    '''

    def __init__(self, message=LverageErrorCode.INVALID_EMAIL.value):
        super().__init__(message, LverageErrorCode.INVALID_EMAIL)


class IdentityThresholdError(LverageError):
    '''
    Exception raised for invalid identity threshold.
    '''

    def __init__(self, message=LverageErrorCode.INVALID_IDENTITY_THRESHOLD.value):
        super().__init__(message, LverageErrorCode.INVALID_IDENTITY_THRESHOLD)

class CantConvertToTaxIDError(LverageError):
    '''
    Exception raised for errors in converting ortholog to taxon ID.
    '''

    def __init__(self, message=LverageErrorCode.CANT_CONVERT_TO_TAX_ID.value):
        super().__init__(message, LverageErrorCode.CANT_CONVERT_TO_TAX_ID)

class BlastHitCountError(LverageError):
    '''
    Exception raised for invalid number of hits for BlastP.
    '''

    def __init__(self, message=LverageErrorCode.INVALID_BLAST_HIT_COUNT.value):
        super().__init__(message, LverageErrorCode.INVALID_BLAST_HIT_COUNT)

class MotifHitCountError(LverageError):
    '''
    Exception raised for invalid number of hits for motif database.
    '''

    def __init__(self, message=LverageErrorCode.INVALID_MOTIF_HIT_COUNT.value):
        super().__init__(message, LverageErrorCode.INVALID_MOTIF_HIT_COUNT)

class EScoreThresholdError(LverageError):
    '''
    Exception raised for invalid e-score threshold.
    '''

    def __init__(self, message=LverageErrorCode.INVALID_ESCORE_THRESHOLD.value):
        super().__init__(message, LverageErrorCode.INVALID_ESCORE_THRESHOLD)