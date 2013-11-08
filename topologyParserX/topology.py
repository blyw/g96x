infile = open('FUL.top','r')

#titleBlock = False
#physicalConstantsBlock = False
#versionBlock = False
#atomTypeNameBlock = False
#residueNameBlock = False
#soluteAtomBlock = False
#bondStretchTypeBlock = False
#bondBlock = False
#bondHBlock = False
#bondAngleBendType = False
#bondAngleHydrogen = False
#bondAngle = False
#improperDihedralTypeBlock = False
#improperDihedralHydrogenBlock = False
#improperDihedralBlock = False
#torsionDihedralTypeBlock = False
#dihedralHydrogenBlock = False
#dihedralBlock = False
#crossDihedralHydrogenBlock = False
#crossDihedralBlock = False
#ljParametersBlock = False
#soluteMoleculesBlock = False
#temperatureGroupsBlock = False
#pressureGroupsBlock = False
#ljExceptionsBlock = False
#solventAtomBlock = False
#solventConstraintsBlock = False
#BONDANGLEBENDTYPE

blockNames = ('TITLE','PHYSICALCONSTANTS','TOPVERSION','ATOMTYPENAME','RESNAME','SOLUTEATOM',
              'BONDSTRETCHTYPE','BONDH','BOND','BONDANGLEBENDTYPE','BONDANGLEH','BONDANGLE',
              'IMPDIHEDRALTYPE','IMPDIHEDRALH','IMPDIHEDRAL','TORSDIHEDRALTYPE','DIHEDRALH',
              'DIHEDRAL','CROSSDIHEDRALH','CROSSDIHEDRAL','LJPARAMETERS','SOLUTEMOLECULES',
              'TEMPERATUREGROUPS','PRESSUREGROUPS','LJEXCEPTIONS','SOLVENTATOM','SOLVENTCONSTR')

data = {}
for blockName in blockNames:
    data[blockName] = []

endOfBlock = True
currentBlockName = ''

for line in infile:
    line = line.strip()
    if len(line) > 0:
        if line[0] != '#':
            if endOfBlock == False and line[:3]=='END':
                endOfBlock = True
            elif endOfBlock == True:
                for blockName in blockNames:
                    if blockName==line:
                       currentBlockName = blockName
                       endOfBlock = False
                       break
                print line
            else:
                data[currentBlockName] += line.split()

outfile = open('FUL.xml','w')
outfile.writelines('<?xml version="1.0" encoding="UTF-8"?>\n<grx>\n <topology version="%s">\n' % tuple(data[blockNames[2]]))
#title block
temp = ''
for item in data[blockNames[0]]:
    temp += item + ' '
outfile.writelines('  <title><![CDATA[%s]]></title>\n' % temp)

#physical constants
outfile.writelines((
                   '  <physical_constants>\n' +
                   '   <permittivity description="FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)">%s</permittivity>\n' +
                   '   <plank description="HBAR: Planck\'s constant HBAR = H/(2* PI)">%s</planck>\n' +
                   '   <light_speed description="SPDL: Speed of light (nm/ps)">%s</light_speed>\n' +
                   '   <boltzman description="BOLTZ: Boltzmann\'s constant kB">%s</boltzman>\n' +
                   '  </physical_constants>\n') % tuple(data[blockNames[1]]))
#atomtypename
outfile.writelines('  <atom_types count="%s">\n' % data[blockNames[3]][0].strip())
i = 1
for item in data[blockNames[3]][1:]:
    outfile.writelines('   <type id="'+str(i)+'">'+ item.strip()+ '</type>\n')
    i+=1
outfile.writelines('  </atom_types>\n')

#resname
outfile.writelines('  <residues count="%s">\n' % data[blockNames[4]][0].strip())
i = 1
for item in data[blockNames[4]][1:]:
    outfile.writelines('   <name id="'+str(i)+'">'+ item.strip()+ '</name>\n')
    i+=1
outfile.writelines('  </residues>\n')

#soluteatom
atom_counter = 1
i = 0
charge_group = 1
outfile.writelines('  <solute_atoms count="%s">\n' % data[blockNames[5]][0])
while i < int(data[blockNames[5]][0]):
    charge_group_indicator = int(data[blockNames[5]][atom_counter+6])
    outfile.writelines('   <atom id="%s" residue_id="%s" name="%s" type="%s" mass="%s" charge="%s" charge_group="%s">\n' % (tuple(data[blockNames[5]][atom_counter:atom_counter+6]) + (charge_group,)))
    atom_counter += 7
    exclusion_counter = int(data[blockNames[5]][atom_counter])
    outfile.writelines('    <exclusions>\n')
    j = 1
    while j <= exclusion_counter:
     atom_counter += 1
     outfile.writelines('     <exclusion>'+data[blockNames[5]][atom_counter]+'</exclusion>\n')
     j+=1
    atom_counter+=1
    exclusion_counter = int(data[blockNames[5]][atom_counter])
    j = 1
    while j <= exclusion_counter:
     atom_counter += 1
     outfile.writelines('     <exclusion_14>'+data[blockNames[5]][atom_counter]+'</exclusion_14>\n')
     j+=1
    atom_counter+=1
    outfile.writelines('    </exclusion>\n')
    outfile.writelines('   </atom>\n')
    i += 1
    if charge_group_indicator==1:
        charge_group += 1
outfile.writelines('  </solute_atoms>\n')

#bond_types
outfile.writelines('  <bond_types count="%s">\n' % data[blockNames[6]][0].strip())
i = 1
while i/3 < int(data[blockNames[6]][0]):
    outfile.writelines('   <type id="%d" quadratic="%s" harmonic="%s" length="%s"/>\n' % (i/3+1, data[blockNames[6]][i], data[blockNames[6]][i+1], data[blockNames[6]][i+2]))
    i+=3
outfile.writelines('  </bond_types>\n')

#bond_hydrogen
outfile.writelines('  <bond_hydrogen count="%s">\n' % data[blockNames[7]][0].strip())
i = 1
while i/3 < int(data[blockNames[7]][0]):
    outfile.writelines('   <type id="%d" atom_1="%s" atom_2="%s" type="%s"/>\n' % (i/3+1, data[blockNames[7]][i], data[blockNames[7]][i+1], data[blockNames[7]][i+2]))
    i+=3
outfile.writelines('  </bond_hydrogen>\n')

#bond_hydrogen
outfile.writelines('  <bond_nonhydrogen count="%s">\n' % data[blockNames[8]][0].strip())
i = 1
while i/3 < int(data[blockNames[8]][0]):
    outfile.writelines('   <type id="%d" atom_1="%s" atom_2="%s" type="%s"/>\n' % (i/3+1, data[blockNames[8]][i], data[blockNames[8]][i+1], data[blockNames[8]][i+2]))
    i+=3
outfile.writelines('  </bond_nonhydrogen>\n')

#bondanglebendtype
outfile.writelines('  <bond_angle_types count="%s">\n' % data[blockNames[9]][0].strip())
i = 1
while i/3 < int(data[blockNames[9]][0]):
    outfile.writelines('   <type id="%d" force_constant_cosine="%s" force_constant="%s" angle="%s"/>\n' % (i/3+1, data[blockNames[9]][i], data[blockNames[9]][i+1], data[blockNames[9]][i+2]))
    i+=3
outfile.writelines('  </bond_angle_types>\n')

#bond_angle_hydrogen
outfile.writelines('  <bond_angle_hydrogen count="%s">\n' % data[blockNames[7]][0].strip())
i = 1
while i/4 < int(data[blockNames[10]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" type="%s"/>\n' % (i/4+1, data[blockNames[10]][i], data[blockNames[10]][i+1], data[blockNames[10]][i+2], data[blockNames[10]][i+3]))
    i+=4
outfile.writelines('  </bond_angle_hydrogen>\n')

#bond_angle_hydrogen
outfile.writelines('  <bond_angle_nonhydrogen count="%s">\n' % data[blockNames[11]][0].strip())
i = 1
while i/4 < int(data[blockNames[11]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" type="%s"/>\n' % (i/4+1, data[blockNames[11]][i], data[blockNames[11]][i+1], data[blockNames[11]][i+2], data[blockNames[11]][i+3]))
    i+=4
outfile.writelines('  </bond_angle_nonhydrogen>\n')

#improper_diheddral_type
outfile.writelines('  <improper_dihedral_types count="%s">\n' % data[blockNames[12]][0].strip())
i = 1
while i/2 < int(data[blockNames[12]][0]):
    outfile.writelines('   <angle id="%d" force_constant="%s" angle="%s"/>\n' % (i/2+1, data[blockNames[12]][i], data[blockNames[12]][i+1]))
    i+=2
outfile.writelines('  </improper_dihedral_types>\n')

#improper_dihedral_hydrogen
outfile.writelines('  <improper_dihedral_hydrogen count="%s">\n' % data[blockNames[13]][0].strip())
i = 1
while i/5 < int(data[blockNames[13]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" atom_4="%s" type="%s"/>\n' % (i/4+1, data[blockNames[13]][i], data[blockNames[13]][i+1], data[blockNames[13]][i+2], data[blockNames[13]][i+3], data[blockNames[13]][i+4]))
    i+=5
outfile.writelines('  </improper_dihedral_hydrogen>\n')

#improper_dihedral_nonhydrogen
outfile.writelines('  <improper_dihedral_nonhydrogen count="%s">\n' % data[blockNames[14]][0].strip())
i = 1
while i/5 < int(data[blockNames[14]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" atom_4="%s" type="%s"/>\n' % (i/4+1, data[blockNames[14]][i], data[blockNames[14]][i+1], data[blockNames[14]][i+2], data[blockNames[14]][i+3], data[blockNames[14]][i+4]))
    i+=5
outfile.writelines('  </improper_dihedral_nonhydrogen>\n')

#torsiondihedraltype
outfile.writelines('  <dihedral_types count="%s">\n' % data[blockNames[15]][0].strip())
i = 1
while i/3 < int(data[blockNames[15]][0]):
    outfile.writelines('   <angle id="%d" force_constant="%s" phase_shift_angle="%s" multiplicity="%s" />\n' % (i/3+1, data[blockNames[15]][i], data[blockNames[15]][i+1], data[blockNames[15]][i+2]))
    i+=3
outfile.writelines('  </dihedral_types>\n')

#dihedral_hydrogen
outfile.writelines('  <dihedral_hydrogen count="%s">\n' % data[blockNames[16]][0].strip())
i = 1
while i/5 < int(data[blockNames[16]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" atom_4="%s" type="%s"/>\n' % (i/4+1, data[blockNames[16]][i], data[blockNames[16]][i+1], data[blockNames[16]][i+2], data[blockNames[16]][i+3], data[blockNames[16]][i+4]))
    i+=5
outfile.writelines('  </dihedral_hydrogen>\n')

#dihedral_nonhydrogen
outfile.writelines('  <dihedral_nonhydrogen count="%s">\n' % data[blockNames[17]][0].strip())
i = 1
while i/5 < int(data[blockNames[17]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" atom_4="%s" type="%s"/>\n' % (i/4+1, data[blockNames[17]][i], data[blockNames[17]][i+1], data[blockNames[17]][i+2], data[blockNames[17]][i+3], data[blockNames[17]][i+4]))
    i+=5
outfile.writelines('  </dihedral_nonhydrogen>\n')

#dihedral_hydrogen
outfile.writelines('  <cross_dihedral_hydrogen count="%s">\n' % data[blockNames[18]][0].strip())
i = 1
while i/9 < int(data[blockNames[18]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" atom_4="%s" atom_5="%s" atom_6="%s" atom_7="%s" atom_8="%s" type="%s"/>\n' % (i/9+1, data[blockNames[18]][i], data[blockNames[18]][i+1], data[blockNames[18]][i+2], data[blockNames[18]][i+3], data[blockNames[18]][i+4], data[blockNames[18]][i+5], data[blockNames[18]][i+6], data[blockNames[18]][i+7], data[blockNames[18]][i+8]))
    i+=9
outfile.writelines('  </cross_dihedral_hydrogen>\n')

#dihedral_nonhydrogen
outfile.writelines('  <cross_dihedral_nonhydrogen count="%s">\n' % data[blockNames[19]][0].strip())
i = 1
while i/9 < int(data[blockNames[19]][0]):
    outfile.writelines('   <angle id="%d" atom_1="%s" atom_2="%s" atom_3="%s" atom_4="%s" atom_5="%s" atom_6="%s" atom_7="%s" atom_8="%s" type="%s"/>\n' % (i/9+1, data[blockNames[19]][i], data[blockNames[19]][i+1], data[blockNames[19]][i+2], data[blockNames[19]][i+3], data[blockNames[19]][i+4], data[blockNames[19]][i+5], data[blockNames[19]][i+6], data[blockNames[19]][i+7], data[blockNames[19]][i+8]))
    i+=9
outfile.writelines('  </cross_dihedral_nonhydrogen>\n')

#dihedral_nonhydrogen
outfile.writelines('  <lj_parameters count="%s">\n' % data[blockNames[20]][0].strip())
i = 1
while i/6 < int(data[blockNames[20]][0]):
    outfile.writelines('   <angle id="%d" atom_type="%s" atom_vwd_type="%s" c12="%s" c6="%s" c12_14="%s" c6_14="%s"/>\n' % (i/6+1, data[blockNames[20]][i], data[blockNames[20]][i+1], data[blockNames[20]][i+2], data[blockNames[20]][i+3], data[blockNames[20]][i+4], data[blockNames[20]][i+5]))
    i+=6
outfile.writelines('  </lj_parameters>\n')

#solutemolecules
outfile.writelines('  <solute_molecules count="%s">\n' % data[blockNames[21]][0].strip())
j = int(data[blockNames[21]][0])
i = 1
while i < int(data[blockNames[21]][0])+1:
    outfile.writelines('   <molecule id="%d" last_atom="%s" />\n' % (i, data[blockNames[21]][i]))
    i+=1
outfile.writelines('  </solute_molecules>\n')

#temperaturegroups
outfile.writelines('  <temperature_groups count="%s">\n' % data[blockNames[22]][0].strip())
j = int(data[blockNames[22]][0])
i = 1
while i < int(data[blockNames[22]][0])+1:
    outfile.writelines('   <group id="%d" last_atom="%s" />\n' % (i, data[blockNames[22]][i]))
    i+=1
outfile.writelines('  </temperature_groups>\n')

#pressuregroups
outfile.writelines('  <pressure_group count="%s">\n' % data[blockNames[23]][0].strip())
j = int(data[blockNames[23]][0])
i = 1
while i < int(data[blockNames[23]][0])+1:
    outfile.writelines('   <group id="%d" last_atom="%s" />\n' % (i, data[blockNames[23]][i]))
    i+=1
outfile.writelines('  </pressure_group>\n')

#ljexception
outfile.writelines('  <lj_exceptions count="%s">\n' % data[blockNames[24]][0].strip())
i = 1
while i/4 < int(data[blockNames[24]][0]):
    outfile.writelines('   <exception id="%d" type1="%s" type2="%s" c12="%s" c6="%s"/>\n' % (i/4+1, data[blockNames[24]][i], data[blockNames[24]][i+1], data[blockNames[24]][i+2], data[blockNames[24]][i+3]))
    i+=4
outfile.writelines('  </lj_exceptions>\n')

#SOLVENTATOM
outfile.writelines('  <solvent_atoms count="%s">\n' % data[blockNames[25]][0].strip())
i = 1
while i/5 < int(data[blockNames[25]][0]):
    outfile.writelines('   <atom id="%d" number="%s" name="%s" type="%s" mass="%s" charge="%s"/>\n' % (i/5+1, data[blockNames[25]][i], data[blockNames[25]][i+1], data[blockNames[25]][i+2], data[blockNames[25]][i+3], data[blockNames[25]][i+4]))
    i+=5
outfile.writelines('  </solvent_atoms>\n')

#SOLVENTCONSTR
outfile.writelines('  <solvent_constraints count="%s">\n' % data[blockNames[26]][0].strip())
i = 1
while i/3 < int(data[blockNames[26]][0]):
    outfile.writelines('   <exception id="%d" atom1="%s" atom2="%s" length="%s"/>\n' % (i/3+1, data[blockNames[26]][i], data[blockNames[26]][i+1], data[blockNames[26]][i+2]))
    i+=3
outfile.writelines('  </solvent_constraints>\n')

outfile.writelines(' </topology>\n</grx>\n')