"""
This file contains the template for all OrfSearcher classes. 
Any OrfSearcher tool that is utilized with Lverage MUST inherit from this class!

Copyright (C) <RELEASE_YEAR_HERE> Bradham Lab

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

Correspondence: 
    - Cynthia A. Bradham - cbradham@bu.edu - *
    - Anthony B. Garza   - abgarza@bu.edu  - **
    - Stephanie P. Hao   - sphao@bu.edu    - **
    - Yeting Li          - yetingli@bu.edu - **

    \* Principle Investigator, ** Software Developers
"""

class OrfSearcherTemplate:
    """
    Abstract class for OrfSearcher classes, containing all necessary functionalities
    This should not be instantiated!
    """

    def __init__(self):
        """Constructor"""
        pass

    def search(self, **kwargs):
        """
        Method for searching a DNA sequence for open-reading-frames
        This method should be overridden by a concrete class

        Returns
        -------
        list
            List of str objects representing the open-reading-frames found in the sequence
        
        Raises
        ------
        NotImplementedError
        """

        raise NotImplementedError