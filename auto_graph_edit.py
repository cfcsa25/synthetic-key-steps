from rdkit import Chem
from rdkit.Chem import rdChemReactions
import pprint

def get_stereo_edits(atom, chiral_flags):
    """Record stereo edits for chiral atoms."""
    if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
        atom_map = atom.GetAtomMapNum()
        chiral_flags[atom_map] = 1
        return (atom_map, atom_map, 1)
    return None

def record_bond_change(bond, direction):
    """Return bond change edit based on bond type and direction (+/-)."""
    bond_type = bond[1]
    atom1, atom2 = bond[0]
    bond_order_map = {
        Chem.rdchem.BondType.SINGLE: 1,
        Chem.rdchem.BondType.DOUBLE: 2,
        Chem.rdchem.BondType.TRIPLE: 3,
        Chem.rdchem.BondType.AROMATIC: 1.5,
        Chem.rdchem.BondType.QUADRUPLE: 4
    }
    if bond_type not in bond_order_map:
        raise ValueError("Warning: bond type not recorded.")
    return (atom1, atom2, direction * bond_order_map[bond_type])



def _bonds_from_template(template_mol, kekulize=False):
    """
    Extract bonds between mapped atoms from a template molecule.
    Returns: dict {(atom1_mapNum, atom2_mapNum): rdkit.Chem.rdchem.BondType}
    """
    # work on a copy to optionally kekulize safely

    copy_mol = Chem.Mol(template_mol)

    if kekulize:
        try:
            Chem.Kekulize(copy_mol, clearAromaticFlags=True)
        except Exception:
            # fallback if kekulization fails; just keep as is
            pass

    idx_to_map = {}
    for a in copy_mol.GetAtoms():
        m = a.GetAtomMapNum()
        if m:  # >0 means mapped, since we don't really care the unmapped atoms (neglected in matrix)
            idx_to_map[a.GetIdx()] = m

    edge2type = {}
    for b in copy_mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a1 in idx_to_map and a2 in idx_to_map:
            i, j = sorted((idx_to_map[a1], idx_to_map[a2])) # record the bond mapping number instead of the atom index
            edge2type[(i, j)] = b.GetBondType()
    return edge2type

def summarize_bond_changes(rxn, kekulize=False):
    """
    Given an rdChemReactions.ChemicalReaction, report which bonds are formed,
    broken, or change order, based on template differences.
    """
    # union across all reactant templates
    r_edge2type = {}
    for i in range(rxn.GetNumReactantTemplates()):
        tmpl = rxn.GetReactantTemplate(i) #rdkit.Chem.rdchem.Mol object
        r_edge2type.update(_bonds_from_template(tmpl, kekulize=kekulize))

    # union across all product templates
    p_edge2type = {}
    # print(rxn.GetNumProductTemplates())
    for i in range(rxn.GetNumProductTemplates()):
        tmpl = rxn.GetProductTemplate(i) #rdkit.Chem.rdchem.Mol object
        p_edge2type.update(_bonds_from_template(tmpl, kekulize=kekulize))

    # print(r_edge2type)

    r_edges = set(r_edge2type)
    p_edges = set(p_edge2type)

    # print(r_edges)
    # print(p_edges)

    broken = sorted(r_edges - p_edges)              # present only in reactants
    formed = sorted(p_edges - r_edges)              # present only in products
    common = r_edges & p_edges


    order_changes = [
        (e, r_edge2type[e], p_edge2type[e])
        for e in sorted(common)
        if r_edge2type[e] != p_edge2type[e]
    ]

    return {
        "broken": [(e, r_edge2type[e]) for e in broken],
        "formed": [(e, p_edge2type[e]) for e in formed],
        "order_changes": order_changes,
    }

def process_reaction(rxn_smarts, chiral_flags, pad_info, product_total_atoms):
    """Process a single reaction SMARTS and return its changelog."""
    rxn_obj = rdChemReactions.ReactionFromSmarts(rxn_smarts)
    rxn = rdChemReactions.ChemicalReaction(rxn_obj)
    bond_changes = summarize_bond_changes(rxn_obj)

    log = {'pad': None, 'pad_elem': [], 'edits': []}
    pad_atoms = set()

    if len(rxn.GetProducts()) == 1:
        product = rxn.GetProducts()[0]
        for atom in product.GetAtoms():
            if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                atom_map = atom.GetAtomMapNum()
                if chiral_flags.get(atom_map, 0) != 0:
                    log['edits'].append((atom_map, atom_map, -1))
                    chiral_flags[atom_map] -= 1
    for bond in bond_changes.get("formed", []):
        pad_atoms.update(bond[0])
        log['edits'].append(record_bond_change(bond, -1))

    for bond in bond_changes.get("broken", []):
        pad_atoms.update(bond[0])
        log['edits'].append(record_bond_change(bond, 1))

    for atom_map in pad_atoms:
        if atom_map > product_total_atoms:
            log['pad_elem'].append(pad_info[atom_map])

    log['pad'] = len(log['pad_elem'])
    return log

# === Main Execution ===
changelogs = []
chiral_flags = {}

# Step 1: Initial stereo edits from final product
first_log = {'pad': 0, 'pad_elem': [], 'edits': []}
for atom in target.GetAtoms():
    edit = get_stereo_edits(atom, chiral_flags)
    if edit:
        first_log['edits'].append(edit)

# Step 2: Process reactions in reverse
for rxn_smarts in reversed(rxns_remapped):
    log = process_reaction(rxn_smarts, chiral_flags, pad_info, product_total_atoms)
    changelogs.append(log)

# Final check
if any(value != 0 for value in chiral_flags.values()):
    raise ValueError("Warning: One or more flags are not set to 0.")

# Final changelog list
changelogs = [first_log] + changelogs[::-1]
pprint.pprint(changelogs)
