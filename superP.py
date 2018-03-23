#!/bin/env python
from Bio.PDB import *
import numpy as np
import argparse


def superpose(chain1,res1,chain2,res2,structure):
    '''Superposes the atoms in chain2 res2 over those in chain1 res1
    Returns a tuple containing the rotation matrix and translation vector
    needed to apply to chain2 res2'''
    # Perform superposition
    sup = Superimposer()
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
    # print([atom for atom in atoms])
    i = 0
    size = len(atoms)
    while i < size:
        if atoms[i].get_id() not in phosphate:
            # print("Popping",atoms[i].get_id())
            atoms.pop(i)
            size -= 1
            # print([atom for atom in atoms])
        else:
            i += 1
    # Make sure they're all in the right order
    atoms = sorted(atoms,key=lambda x: phosphate.index(x.get_id()))
    return atoms
    
def readMatrix(matfile):
    min = open(matfile)
    output = [line.split(',') for line in min]
    min.close()
    output.pop(0)
    return output
    
def removeSimilar(pdb,mat1,mat2,cutoff):
    # initialize output
    output = []
    # for every entry in mat1
    for entry1 in mat1:
        # for every entry in mat2
        failed = False
        for entry2 in mat2:
            print(entry1[0:4],entry2[0:4],'        ',end='\r')
            # load pdb
            parser = PDBParser()
            struc = parser.get_structure(pdb,pdb)
            # extract mat1 atoms
            # print(entry1[0],entry1[1])
            # print([atom for atom in struc[0][entry1[0]][int(entry1[1])]])
            atoms1 = extractPhosphate([atom for atom in struc[0][entry1[0]][int(entry1[1])]])
            # perform mat2 sup on mat1 atoms
            rot = [float(x) for x in entry2[4:13]]
            rot = np.reshape(np.array(rot),(3,3))
            tran = [float(x) for x in entry2[13:16]]
            tran = np.array(tran)
            # if dist(mat1 a,b) < cutoff, remove from output
            # print([atom.get_id() for atom in atoms1])
            atoms1P = [atom for atom in atoms1 if atom.get_id()=='P']
            atoms2 =[atom for atom in struc[0][entry1[2]][int(entry1[3])] if atom.get_id()=='P']
            atoms1P[0].transform(rot,tran)
            dist = atoms2[0] - atoms1P[0]
            if dist < cutoff:
                failed = True
                print("\nRemoving",entry1[0:4])
                break
        if not failed:
            output.append(entry1)
    # return output
    return output
    
def main():
    parser = argparse.ArgumentParser()
    # input pdb
    parser.add_argument('-f',dest='pdb',help='DNA PDB file to be superposed')
    # output csv file
    parser.add_argument('-o',dest='outfile',help='CSV file to save output to')
    args = parser.parse_args()
    pdb = args.pdb
    outfile = args.outfile
    
    # Read in pdb file
    pdbParser = PDBParser()
    structure = pdbParser.get_structure('structure',pdb)
    # Initialize output list
    output = []
    # Iterating over every chain and residue
    for chain1 in structure[0]:
        for res1 in chain1:
            for chain2 in structure[0]:
                for res2 in chain2:
                    print(chain1.get_id(),res1.get_id()[1],chain2.get_id(),res2.get_id()[1],end="\r")
                    # Superpose their phosphates and add them to the list
                    sup = superpose(chain1,res1,chain2,res2,structure)
                    if sup != -1:
                        output.append((chain1.get_id(),res1.get_id()[1],chain2.get_id(),res2.get_id()[1],sup))
    # Write out the csv file
    writeCSV(output,outfile)
    
if __name__ == "__main__": main()