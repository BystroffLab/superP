#!/bin/env python
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Model
import optparse


def checkFirst(res1,chain1,res2,chain2):
    if res1 == chain1.get_list()[0] and res2 != chain2.get_list()[0]: return True
    return False

def superpose(chain1,res1,chain2,res2,structure):
    
    #Perform superposition
    sup = Bio.PDB.Superimposer()
    fixed  = [atom for atom in res1]
    moving = [atom for atom in res2]
    
    #Get rid of hydrogens which may or may not exist in either structure and keep only
    #the phosphate backbone
    extractBackbone(fixed)
    extractBackbone(moving)
    
    #If either starts with the first nucleotide, there will be no phosphate at the beginning
    #Fine if both do!
    if checkFirst(res2,chain2,res1,chain1):
        for i in range(3):
            fixed.pop(0)
    if checkFirst(res1,chain1,res2,chain2):
        for i in range(3):
            moving.pop(0)
    #Calculate the superimposition matrices
    try:
        sup.set_atoms(fixed,moving)
    except Exception as e:
        print ("Error superimposing PX.")
        print ('fixed atoms:',len(fixed),fixed)
        print ('moving atoms:',len(moving),moving)
        raise e
    #return rotation matrix and translation vector
    return sup.rotran
    
def writeCSV(matrices,outfile):
    '''
    chain1,res1,chain2,res2,m11,m12,m13,m21,m22,m23,m31,m32,m33,v1,v2,v3
    '''
    out = open(outfile,"w+")
    out.write("chain1,res1,chain2,res2,m11,m12,m13,m21,m22,m23,m31,m32,m33,v1,v2,v3\n")
    for (chain1,res1,chain2,res2,(rot,tran)) in matrices:
        out.write("%s,%i,%s,%i,"%(chain1,res1,chain2,res2))
        out.write("%f,%f,%f,"%(rot[0][0],rot[0][1],rot[0][2]))
        out.write("%f,%f,%f,"%(rot[1][0],rot[1][1],rot[1][2]))
        out.write("%f,%f,%f,"%(rot[2][0],rot[2][1],rot[2][2]))
        out.write("%f,%f,%f\n"%(tran[0],tran[1],tran[2]))
    out.close()
        
    
def extractBackbone(atoms):
    #DNA backbone
    backbone = ['P','OP1','OP2','O5\'','C5\'','C4\'','O4\'','C3\'','C4\'','O4\'','C2\'','C1\'']
    i = 0
    size = len(atoms)
    while i < size:
        if atoms[i].get_id() not in backbone:
            atoms.pop(i)
            size -= 1
        else:
            i += 1
            
def main():
    parser = optparse.OptionParser()
    #input pdb
    parser.add_option('-f',dest='pdb',help='PDB file to be superposed')
    #output csv file
    parser.add_option('-o',dest='outfile',help='CSV file to save output to')
    
    (options,args) = parser.parse_args()
    pdb = options.pdb
    outfile = options.outfile
    
    pdbParser = Bio.PDB.PDBParser()
    structure = pdbParser.get_structure('structure',pdb)
    output = []
    for chain1 in structure[0]:
        for res1 in chain1:
            for chain2 in structure[0]:
                for res2 in chain2:
                    print(chain1.get_id(),res1.get_id()[1],chain2.get_id(),res2.get_id()[1])
                    output.append((chain1.get_id(),res1.get_id()[1],chain2.get_id(),res2.get_id()[1],superpose(chain1,res1,chain2,res2,structure)))
    writeCSV(output,outfile)
    
if __name__ == "__main__": main()