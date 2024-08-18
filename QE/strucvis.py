from ase import Atoms
from ase.io import read, write
from tempfile import NamedTemporaryFile
from IPython.display import HTML

class qe2xyz:
    def __init__(self, input_file):
        self.input_file = input_file

    def convert(self, output_file):
        atoms = read(self.input_file, format='espresso-in')
        write(output_file, atoms, format='xyz')

class StructureVisualizer:
    def __init__(self, file_path):
        self.file_path = file_path

    def Structure(self):
        'Return the HTML representation of the atoms object as a string'
        with NamedTemporaryFile('r+', suffix='.html') as ntf:
            atoms = read(self.file_path)
            atoms.center()  # Center the atoms
            atoms.write(ntf.name, format='html')
            ntf.seek(0)
            html = ntf.read()
        return HTML(html)  # Return the HTML string
