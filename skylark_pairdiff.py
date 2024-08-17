import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import joblib

#将SMILES标准化,输入smiles,返回标准化后的smiles,失败返回None
def standardize_smiles(smiles):
    '''
    将SMILES标准化,输入smiles,返回标准化后的smiles,失败返回None
    例:standardize_smiles('C'),返回'C'
    '''
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol)
    except:
        return 'None'

#根据SMILES计算分子量,输入smiles,返回分子量,失败返回0
def get_molecular_weight(smiles):
    '''
    计算分子量,输入smiles,返回分子量,失败返回0
    例:get_molecular_weight('C'),返回12.01
    '''
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Descriptors.MolWt(mol)
    except:
        return 0

#官能团质量变化分组字典，key为质量范围，value为官能团质量变化分组[[被取代官能团，取代官能团],...]
functional_group_mass_dict = joblib.load(r'dict\\model_label_id_dict.joblib')
key_list = list(functional_group_mass_dict.keys())

#获取母体和子体的官能团取代关系，输入为母体和子体的smiles，输出为[被取代官能团，取代官能团]，失败返回'None'
def get_functional_group_change(mother_smiles, son_smiles):
    '''
    判断母体和子体官能团的取代关系,输入母体和子体的smiles,返回[被取代的官能团，取代后的官能团],失败返回'None'
    例:get_functional_group_change('CCCCCC', 'CCCCCCl'),返回['C', 'Cl']
    '''
    try:
        mother_smiles = standardize_smiles(mother_smiles)
        son_smiles = standardize_smiles(son_smiles)
        #若SMILES标准化失败，返回None
        if mother_smiles == 'None' and son_smiles == 'None':
            return 'None'
        else:
            mother_mass = get_molecular_weight(mother_smiles)
            son_mass = get_molecular_weight(son_smiles)
            massdiff = abs(mother_mass - son_mass)
            functional_pair_list = []
            for mass_range in key_list:
                mass_low, mass_high = mass_range.split('_')
                mass_low, mass_high = float(mass_low), float(mass_high)
                if mass_low <= massdiff <= mass_high:
                    functional_pair_dict = functional_group_mass_dict[mass_range]
                    value_list = list(functional_pair_dict.values())
                    for value in value_list:
                        functional1, functional2 = value.split('_')
                        if functional1 != 'unknown':
                            functional_pair_list.append([functional1, functional2])
                    break
            if len(functional_pair_list) == 0:
                return 'None'
            else:
                mother_mol = Chem.MolFromSmiles(mother_smiles)
                son_mol = Chem.MolFromSmiles(son_smiles)
                functional_group_change = []
                for functional_pair in functional_pair_list:
                    replaced_functional_smiles = functional_pair[0]
                    replaced_functional_mol = Chem.MolFromSmiles(replaced_functional_smiles)
                    replace_functional_smiles = functional_pair[1]
                    replace_functional_mol = Chem.MolFromSmiles(replace_functional_smiles)
                    flag1 = mother_mol.HasSubstructMatch(replaced_functional_mol, useChirality=True)
                    flag2 = son_mol.HasSubstructMatch(replace_functional_mol, useChirality=True)
                    if flag1 and flag2:
                        possible_mol_list = AllChem.ReplaceSubstructs(mother_mol, 
                                                                      replaced_functional_mol, 
                                                                      replace_functional_mol)
                        for possible_mol in possible_mol_list:
                            possible_smi = Chem.MolToSmiles(possible_mol)
                            if possible_smi == son_smiles:
                                functional_group_change.append(functional_pair)
                if len(functional_group_change) == 0:
                    return 'None'
                else:
                     #选取字符串长度最短的候选结构
                    min_len = len(functional_group_change[0][0])+len(functional_group_change[0][1])
                    min_len_index = 0
                    for i in range(len(functional_group_change)):
                        if len(functional_group_change[i][0])+len(functional_group_change[i][1]) < min_len:
                            min_len = len(functional_group_change[i][0])+len(functional_group_change[i][1])
                            min_len_index = i
                    return functional_group_change[min_len_index]
    except:
        return 'None'