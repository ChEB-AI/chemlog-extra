from chemlog.base_classifier import Classifier
from chemlog.preprocessing.chebi_data import ChEBIData
from rdkit import Chem

class XMolecularEntityClassifier(Classifier):

    def __init__(self, chebi_version: int = 241):
        super().__init__()
        self.chebi_version = chebi_version
        self.element_class_mapping = self.build_class_element_mapping(chebi_version=chebi_version)
        print(f"Found {len(self.element_class_mapping)} ... molecular entity classes")

    def classify(self, mol_list):
        res = []
        for mol in mol_list:
            res.append([self.element_class_mapping[element_num] for element_num in
                        list(set([atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 0]))
                        if element_num in self.element_class_mapping])
        return res

    def build_class_element_mapping(self, chebi_version: int = 241):
        element_name_to_num = {Chem.GetPeriodicTable().GetElementName(i).lower(): i for i in range(1, 119)}
        data = ChEBIData(chebi_version)
        element_class_mapping = {}
        for chebi_id, properties in data.chebi.items():
            if "name" in properties:
                if " molecular entity" in properties["name"]:
                    element_name = properties["name"].split(" ")[0]
                    if element_name == "organic":
                        element_name = "carbon"
                    if element_name in element_name_to_num:
                        element_class_mapping[element_name_to_num[element_name]] = chebi_id
        return element_class_mapping


class OrganoXCompoundClassifier(Classifier):

    def __init__(self, chebi_version: int = 241):
        super().__init__()
        self.element_class_mapping = self.build_class_element_mapping(chebi_version=chebi_version)
        print(f"Found {len(self.element_class_mapping)} organo-... compound classes")

    def classify(self, mol_list):
        res = []
        for mol in mol_list:
            res.append([self.element_class_mapping[element_num] for element_num in
                        list(set([atom.GetAtomicNum() for atom in mol.GetAtoms()
                                  if any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())]))
                        if element_num in self.element_class_mapping])
        return res

    def build_class_element_mapping(self, chebi_version: int = 241):
        element_name_to_num = {Chem.GetPeriodicTable().GetElementName(i).lower(): i for i in range(1, 119)}
        data = ChEBIData(chebi_version)
        element_class_mapping = {}
        for chebi_id, properties in data.chebi.items():
            if "name" in properties:
                if properties["name"].startswith("organo") and " compound" in properties["name"]:
                    element_name = properties["name"][6:].split(" ")[0]
                    if element_name in element_name_to_num:
                        element_class_mapping[element_name_to_num[element_name]] = chebi_id
        return element_class_mapping

if __name__ == "__main__":
    classifier = XMolecularEntityClassifier()
    mapping = classifier.build_class_element_mapping(chebi_version=241)
    print(classifier.classify([Chem.MolFromSmiles("C12=C(N(C(=O)N(C)C1=O)C)N=CN2C.C1(=CC=CC=C1)C(=O)[O-].[Na+]")]))
    classifier = OrganoXCompoundClassifier()
    print(classifier.classify([Chem.MolFromSmiles("C12=C(N(C(=O)N(C)C1=O)C)N=CN2C.C1(=CC=CC=C1)C(=O)[O-].[Na+]")]))