#!/bin/env python
from __future__ import print_function
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Model
import optparse


def superpose(chain1,res1,chain2,res2,structure):
    '''Superposes the atoms in chain2 res2 over those in chain1 res1
    Returns a tuple containing the rotation matrix and translation vector
    needed to apply to chain2 res2'''
    # Perform superposition
    sup = Bio.PDB.Superimposer()
    fixed  = [atom for atom in res1]
    moving = [atom for atom in res2]
    
    # Get rid of hydrogens which may or may not exist in either structure and keep only
    # the phosphate backbone
    extractPhosphate(fixed)
    extractPhosphate(moving)
    
    # There should be four atoms in each atom set, if not, skip.
    if len(fixed) != 4 or len(moving) != 4: return -1
    # Calculate the rotation matrix and translation vector
    try:
        sup.set_atoms(fixed,moving)
    except Exception as e:
        # Is it due to size mismatch?
        print ("Error superimposing PX.")
        print ('fixed atoms:',len(fixed),fixed)
        print ('moving atoms:',len(moving),moving)
        raise e
    # Return rotation matrix and translation vector
    return (sup.rotran,sup.rms)
    
def writeCSV(matrices,outfile):
    '''Takes our rotation matrices and translation vectors and outputs them as
     a csv file containing the chains,residues, and unrolled matrices and
      vectors in the following format:
    chain1,res1,chain2,res2,m11,m12,m13,m21,m22,m23,m31,m32,m33,v1,v2,v3,rms
    '''
    out = open(outfile,"w+")
    out.write("chain1,res1,chain2,res2,m11,m12,m13,m21,m22,m23,m31,m32,m33,v1,v2,v3,rms\n")
    for (chain1,res1,chain2,res2,((rot,tran),rms)) in matrices:
        out.write("%s,%i,%s,%i,"%(chain1,res1,chain2,res2))
        out.write("%f,%f,%f,"%(rot[0][0],rot[0][1],rot[0][2]))
        out.write("%f,%f,%f,"%(rot[1][0],rot[1][1],rot[1][2]))
        out.write("%f,%f,%f,"%(rot[2][0],rot[2][1],rot[2][2]))
        out.write("%f,%f,%f,%f\n"%(tran[0],tran[1],tran[2],rms))
    out.close()

def extractPhosphate(atoms):
    '''Given an input set of atoms, removes non-phosphate atoms'''
    # DNA backbone
    # backbone = ['P','OP1','OP2','O5\'','C5\'','C4\'','O4\'','C3\'','C4\'','O4\'','C2\'','C1\'']
    phosphate = ['P','OP1','OP2','O5\'']
    i = 0
    size = len(atoms)
    while i < size:
        if atoms[i].get_id() not in phosphate:
            atoms.pop(i)
            size -= 1
        else:
            i += 1
    # Make sure they're all in the right order
    atoms = sorted(atoms,key=lambda x: phosphate.index(x.get_id()))
            
def main():
    parser = optparse.OptionParser()
    # input pdb
    parser.add_option('-f',dest='pdb',help='DNA PDB file to be superposed')
    # output csv file
    parser.add_option('-o',dest='outfile',help='CSV file to save output to')
    (options,args) = parser.parse_args()
    pdb = options.pdb
    outfile = options.outfile
    
    # Read in pdb file
    pdbParser = Bio.PDB.PDBParser()
    structure = pdbParser.get_structure('structure',pdb)
    # Initialize output list
    output = []
    # Iterating over every chain and residue
    for chain1 in structure[0]:
        for res1 in chain1:
            for chain2 in structure[0]:
                for res2 in chain2:
                    print(chain1.get_id(),res1.get_id()[1],chain2.get_id(),res2.get_id()[1])
                    # Superpose their phosphates and add them to the list
                    sup = superpose(chain1,res1,chain2,res2,structure)
                    if sup != -1:
                        output.append((chain1.get_id(),res1.get_id()[1],chain2.get_id(),res2.get_id()[1],sup))
    # Write out the csv file
    writeCSV(output,outfile)
    
if __name__ == "__main__": main()