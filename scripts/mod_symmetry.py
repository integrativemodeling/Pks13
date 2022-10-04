#####################################
# Modeling of Pks13 dimer
# using two-fold symmetry
#
# iecheverria - Sali Lab
######################################

import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.mmcif
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.em
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter

import math
import numpy as np
import sys
import os

top_dir='../'

# Create System and State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()
mols = []

all_atomic_0 = []
seqs = IMP.pmi.topology.Sequences(f'{top_dir}/data/MSMEG_Pks13_seq.fasta')
mol = st.create_molecule('Pks13',
                         sequence=seqs['Pks13'],
                         chain_id='A')
atomic = mol.add_structure(f'{top_dir}/data/pks13_KS_apply_ncs_12-coot-6_real_space.centered.pdb',
                           chain_id='A',
                           res_range=[89,1074],
                           offset=0)
all_atomic_0 += atomic
mol.add_representation(atomic,  
                       resolutions=[1],
                       color='red')


# Add domains from AF
doms = [[1,76],[1078,1173],[1238,1356],[1461,1535],[1539,1816]]
colors = ['yellow','dark cyan','dark green','plum','dark magenta']
for k,d in enumerate(doms):
    atomic_dom = mol.add_structure(f'{top_dir}/data/AF2_pks13_monomer.pdb',
                                    chain_id='A',
                                    res_range=d,
                                    offset=0)
    mol.add_representation(atomic_dom,
                       resolutions=[1],
                       color=colors[k])
    all_atomic_0 += atomic_dom

mol.add_representation(mol[:]-all_atomic_0,  
                       resolutions=[10],
                       color='gray')


mols.append(mol)

##################
# Create clone
#################

all_atomic_1 = []
clone = mol.create_copy(chain_id = 'B')
atomic = clone.add_structure(f'{top_dir}/data/pks13_KS_apply_ncs_12-coot-6_real_space.centered.pdb',
                            res_range=[89,1074],
                            chain_id = 'A')
all_atomic_1 += atomic
clone.add_representation(atomic,
                        resolutions = [1],
                        color='blue')
for k,d in enumerate(doms):
    atomic_dom = clone.add_structure(f'{top_dir}/data/AF2_pks13_monomer.pdb',
                                    chain_id='A',
                                    res_range=d,
                                    offset=0)
    clone.add_representation(atomic_dom,
                       resolutions=[1],
                       color=colors[k])
    all_atomic_1 += atomic_dom

if len(clone[:]-all_atomic_1) > 0:                         
    print('Adding missing residues')
    clone.add_representation(clone[:]-all_atomic_1,
                           resolutions=[10],
                            color = 'light gray')

mols.append(clone)



##############################
# Generate mmcif file
##############################
    
if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput()
    po.system.title = ('Integrative structure determination of the Pks13 dimer')
    s.add_protocol_output(po)

##############################
# Build system
##############################

print('Molecules:', mols)

hier = s.build()


print(IMP.atom.show_with_representations(hier))
# write a single-frame RMF 
out = IMP.pmi.output.Output()
out.init_rmf("structure_sym_ini.rmf3", hierarchies=[hier])
out.write_rmf("structure_sym_ini.rmf3")
    
#########################
# Define DOFs
#########################
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
for i, mol in enumerate(mols):
    dof.create_flexible_beads(mol)
    for j, d in enumerate(doms):
        sel = IMP.atom.Selection(hier,
                                 molecule='Pks13',
                                 residue_indexes=range(d[0],d[1]+1),
                                 copy_index=i).get_selected_particles()
    
        rb_movers,rb = dof.create_rigid_body(sel,
                                             name = f'Pks13_{i}_{j}_flex')
    # Cryo-EM structure
    sel = IMP.atom.Selection(hier,
                                 molecule='Pks13',
                                 residue_indexes=list(range(89,530))+list(range(588,1075)),
                                 copy_index=i).get_selected_particles()

    rb_movers,rb = dof.create_rigid_body(sel,
                                         name = f'Pks13_cryo_{i}')
    if i==0:
        cryo_movers = rb_movers
        print(rb_movers)



# Create a symmetry constraint
#  A constrant is invariant: IMP will automatically move all clones to
#  match the reference
#  If instead you want some more flexiblity, consider
#  IMP.pmi.restraints.stereochemistry.SymmetryRestraint
center = IMP.algebra.Vector3D([0, 1, 0])

fit_t = IMP.algebra.Transformation3D(
        IMP.algebra.get_rotation_from_matrix(
            -0.979135,   0.202447, -0.0176007, 
             0.202344,   0.96331,  -0.176324, 
            -0.0187413, -0.176207, -0.984175 ),
        IMP.algebra.Vector3D(0.0167301, 0.00551412, 0.00550137))

dof.constrain_symmetry(mols[0], mols[1], fit_t)
mdl.update()  # propagates coordinates

# write a single-frame RMF to view the helix
out = IMP.pmi.output.Output()
out.init_rmf("test_symmetry.rmf3", hierarchies=[hier])
out.write_rmf("test_symmetry.rmf3")

output_objects = []

#################################
# Connectivity restraint
#################################

cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mols[0])
cr.set_label('Pks13.0')
cr.add_to_model()
output_objects.append(cr)

#################################
# Excluded volume restraint
#################################

evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols[0])
evr1.add_to_model()
evr1.set_weight(1.0)
evr1.set_label('intra')
output_objects.append(evr1)

evr2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols[0],
                                                               other_objects=mols[1])
evr2.add_to_model()
evr2.set_weight(1.0)
evr2.set_label('inter')
output_objects.append(evr2)

print('Excluded volume:', evr1.evaluate(), evr2.evaluate())

#################################
# Crosslinking restraint
#################################
rmf_restraints = []

# INITIALIZE DB    

cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein1")
cldbkc.set_protein2_key("Protein2")
cldbkc.set_residue1_key("AbsPos1")
cldbkc.set_residue2_key("AbsPos2")
cldbkc.set_unique_id_key("Id")
cldbkc.set_psi_key("Score")

# XLs RESTRAINT
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file(f"{top_dir}/data/Interlinks_Pks13_XLEHP1007-15_1064-69_filtered_nonambiguos_modeling_unique_modeling.csv")

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                            database=cldb,
                                                                            resolution=1.0,
                                                                            length=26.0,
                                                                            slope=0.00001)
xl1.add_to_model()
xl1.set_weight(10.0)

rmf_restraints.append(xl1)
output_objects.append(xl1)

print(xl1.get_particles_to_sample())
print('XLs score:', xl1.evaluate())


#################################
# Randomize configurations
#################################

sel = IMP.atom.Selection(hier,
                         molecule='Pks13',
                         residue_indexes=list(range(1,76))+list(range(1075,1816))).get_selected_particles()

IMP.pmi.tools.shuffle_configuration(sel,
                                    bounding_box=((-200, -50, -200), (200, 200, 200)),
                                    avoidcollision_rb=False)

print(len(dof.get_movers()), cryo_movers)
sel_movers = [m for m in dof.get_movers() if m not in cryo_movers]
print(len(sel_movers))

#################################
# Sampling
#################################

num_frames = 30000
if '--mmcif' in sys.argv or '--test' in sys.argv:
    num_frames=5

mc = IMP.pmi.macros.ReplicaExchange0(mdl,
                                      root_hier=hier,                       
                                      crosslink_restraints=rmf_restraints,       
                                      monte_carlo_sample_objects=sel_movers,  
                                      global_output_directory="output/",
                                      output_objects=output_objects,
                                      replica_exchange_maximum_temperature=3.0,
                                      monte_carlo_steps=10,
                                      number_of_frames = num_frames,
                                      number_of_best_scoring_models=0)

mc.execute_macro()



##############################
# Generate mmcif
##############################  
if '--mmcif' in sys.argv:
    import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    
    fname = os.path.join(top_dir,"data/Interlinks_Pks13_XLEHP1007-15_1064-69_filtered_nonambiguos_modeling_unique_modeling.csv")

    # Add publication
    po.system.title = "Structure and dynamics of the essential endogenous mycobacterial polyketide synthase Pks13"
    #po.system.citations.append(ihm.Citation.from_pubmed_id(32719457))
    
    s = po.system
    print("restraint datasets:", [r.dataset for r in s.restraints])
    # Datasets for XL-MS restraint
    for r in s.restraints:
        if isinstance(r, ihm.restraint.CrossLinkRestraint):
            r.linker = ihm.cross_linkers.dsso
            print("XL-MS dataset at:", r.dataset.location.path)
            print("Details:", r.dataset.location.details)
    
    # Correct number of output models to account for multiple runs
    protocol = s.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 2500000

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 9999 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=200000,
                            num_models_end=9999))
    

    # Create an ensemble for the cluster
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0", num_models=9999,
                                drmsd=8.3, num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    #rh = RMF.open_rmf_file_read_only('../results/clustering_rigid/cluster.0/cluster_center_model.rmf3')
    #IMP.rmf.link_hierarchies(rh, [hier])
    #IMP.rmf.load_frame(rh, RMF.FrameID(0))
    #del rh
    #model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        for s in np.arange(0,6,1):
            name = asym.split('.')[0]
            fname = f'../results/clustering/cluster.0/LPD_{name}_0_{s}.mrc'
            print('fname', fname)
            loc = ihm.location.OutputFileLocation(fname)
            den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
            # Add to ensemble
            e.densities.append(den)
    
    # Add uniprot of proteins

    lpep = ihm.LPeptideAlphabet()
    d = 'Construct used for cryo-EM'
    sd_pks13 = [ihm.reference.SeqDif(883, lpep['D'], lpep['A'], details=d),
              ihm.reference.SeqDif(884, lpep['A'], lpep['I'], details=d),
              ihm.reference.SeqDif(886, lpep['I'], lpep['A'], details=d)]

    # name : (uniprot id, mutations, [[db_begin, db_end, entity_begin, entity_end]]
    Uniprot={'Pks13.0': ('A0R617',sd_pks13,[[1,1815,1,1815]]),
             'Pks13.1': ('A0R617',sd_pks13,[[1,1815,1,1815]])}

    for prot, (entry, sd, limits) in Uniprot.items():
        print(prot, entry, sd, limits)
        ref = ihm.reference.UniProtSequence.from_accession(entry)
        for seg in limits:
            ref.alignments.append(ihm.reference.Alignment(
                db_begin=seg[0], db_end=seg[1], entity_begin=seg[2], entity_end=seg[3], seq_dif=sd))
            
        po.asym_units[prot].entity.references.append(ref)

    # Point to the raw mass spec data and peaklists used to derive the crosslinks.
    #l = ihm.location.PRIDELocation('PXD019338',
    #                               details='All raw mass spectrometry files and '
    #                               'peaklists used in the study')
    #xl1.dataset.add_primary(ihm.dataset.MassSpecDataset(location=l))
        

    # Replace local links with DOIs
    repos = []
    #for subdir, zipname in make_archive.ARCHIVES.items():
    #    print('subdir', subdir)
    #    repos.append(ihm.location.Repository(
    #        doi="10.5281/zenodo.3836213", root="../%s" % subdir,
    #        url="https://zenodo.org/record/3836213/files/%s.zip" % zipname,
    #        top_directory=os.path.basename(subdir)))
    
    po.system.update_locations_in_repositories(repos)

    with open('IM_Pks13_dimer.cif', 'w') as fh:
        ihm.dumper.write(fh, [po.system])
    #os.unlink('IM_Pks13_dimer.cif')

exit()

