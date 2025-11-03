import os
import csv
from chemlog.base_classifier import Classifier
from chemlog.preprocessing.chebi_data import ChEBIData
from rdkit import Chem

class ExtraClassifier(Classifier):

    def __init__(self, chebi_graph, chebi_version: int = 244):
        super().__init__()
        self.chebi_graph = chebi_graph
        self.chebi_version = chebi_version
        self.element_class_mapping = self.load_element_class_mapping()

    def classify(self, mol_list):
        res = []
        if not isinstance(mol_list, list):
            mol_list = [mol_list]
        for mol in mol_list:
            if not mol:
                res.append({})
                continue
            res.append({cls: self.get_single_classification(mol, element_num)
                        for element_num, cls in self.element_class_mapping.items()})
        return res
    
    def get_single_classification(self, mol, element_num):
        raise NotImplementedError()

    def load_element_class_mapping(self):
        data_path = os.path.join("data", f"chebi_v{self.chebi_version}", f"{self.__class__.__name__}_element_class_mapping.csv")
        if os.path.exists(data_path):
            res = dict()
            with open(data_path, "r") as f:
                reader = csv.reader(f)
                next(reader)  # skip header
                for row in reader:
                    if len(row) == 2:
                        element_num = int(row[0])
                        chebi_id = row[1]
                        res[element_num] = chebi_id
            return res
        else:
            mapping = self.build_class_element_mapping()
            os.makedirs(os.path.dirname(data_path), exist_ok=True)
            with open(data_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["element_num", "chebi_id"])
                for element_num, chebi_id in mapping.items():
                    writer.writerow([element_num, chebi_id])
            return mapping


class XMolecularEntityClassifier(ExtraClassifier):

    
    def build_class_element_mapping(self):
        element_name_to_num = {Chem.GetPeriodicTable().GetElementName(i).lower(): i for i in range(1, 119)}
        element_class_mapping = {}
        for chebi_id, properties in self.chebi_graph.nodes.items():
            if "name" in properties:
                if " molecular entity" in properties["name"]:
                    element_name = properties["name"].split(" ")[0]
                    if element_name == "organic":
                        element_name = "carbon"
                    if element_name in element_name_to_num:
                        element_class_mapping[element_name_to_num[element_name]] = str(chebi_id)
        return element_class_mapping

    def get_single_classification(self, mol, element_num):
        return element_num in list(set([atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 0]))


class OrganoXCompoundClassifier(ExtraClassifier):
        
    def get_single_classification(self, mol, element_num):
        return element_num in list(set([atom.GetAtomicNum() for atom in mol.GetAtoms()
                                  if any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())]))

    def build_class_element_mapping(self):
        element_name_to_num = {Chem.GetPeriodicTable().GetElementName(i).lower(): i for i in range(1, 119)}
        element_class_mapping = {}
        for chebi_id, properties in self.chebi_graph.nodes.items():
            if "name" in properties:
                if properties["name"].startswith("organo") and " compound" in properties["name"]:
                    element_name = properties["name"][6:].split(" ")[0]
                    if element_name in element_name_to_num:
                        element_class_mapping[element_name_to_num[element_name]] = str(chebi_id)
        return element_class_mapping

if __name__ == "__main__":
    from chebifier.utils import build_chebi_graph
    chebi_graph = build_chebi_graph(chebi_version=244)
    classifier = XMolecularEntityClassifier(chebi_graph, chebi_version=244)
    print(classifier.classify([Chem.MolFromSmiles("C12=C(N(C(=O)N(C)C1=O)C)N=CN2C.C1(=CC=CC=C1)C(=O)[O-].[Na+]")]))
    classifier = OrganoXCompoundClassifier(chebi_graph, chebi_version=244)
    print(classifier.classify([Chem.MolFromSmiles("C12=C(N(C(=O)N(C)C1=O)C)N=CN2C.C1(=CC=CC=C1)C(=O)[O-].[Na+]")]))