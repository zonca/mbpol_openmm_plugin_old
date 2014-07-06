import mbpolplugin
from simtk.openmm.app import forcefield
from simtk import unit

#class MBPolElectrostaticsForceGenerator:
#
#    def __init__(self, forceField):
#
#        self.forceField = forceField
#
#    @staticmethod
#    def parseElement(element, forceField):
#
#
#        # <MBPolElectrostaticsForce>
#        #     <Residue class1="OW" class2="HW" class3="HW" thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole-intermolecules="0.055" thole-dipole-dipole-1-2="0.055" thole-dipole-dipole-1-3="0.626" thole-dipole-dipole-2-3="0.626" /> 
#        #     <Atom type="MBPol-O" charge="-5.1966000e-01" damping-factor="0.00131" polarizability="0.00131" />
#        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
#        #     <Atom type="MBPol-M" charge="0" damping-factor="0.00131" polarizability="0" />
#        # </MBPolElectrostaticsForce>
#
#        #   <AmoebaMultipoleForce  direct11Scale="0.0"  direct12Scale="1.0"  direct13Scale="1.0"  direct14Scale="1.0"  mpole12Scale="0.0"  mpole13Scale="0.0"  mpole14Scale="0.4"  mpole15Scale="0.8"  mutual11Scale="1.0"  mutual12Scale="1.0"  mutual13Scale="1.0"  mutual14Scale="1.0"  polar12Scale="0.0"  polar13Scale="0.0"  polar14Intra="0.5"  polar14Scale="1.0"  polar15Scale="1.0"  >
#        # <Multipole class="1"    kz="2"    kx="4"    c0="-0.22620" d1="0.08214" d2="0.00000" d3="0.34883" q11="0.11775" q21="0.00000" q22="-1.02185" q31="-0.17555" q32="0.00000" q33="0.90410"  />
#        # <Multipole class="2"    kz="1"    kx="3"    c0="-0.15245" d1="0.19517" d2="0.00000" d3="0.19687" q11="-0.20677" q21="0.00000" q22="-0.48084" q31="-0.01672" q32="0.00000" q33="0.68761"  />
#
#        generator = MBPolElectrostaticsForceGenerator(forceField)
#
#        forceField._forces.append(generator)
#
#        for residue in element.findall('Residue'):
#            types = forceField._findAtomTypes(atom, 1)
#            if None not in types:
#
#                # k-indices not provided default to 0
#
#                kIndices = [int(atom.attrib['type'])]
#
#                kStrings = [ 'kz', 'kx', 'ky' ]
#                for kString in kStrings:
#                    try:
#                        if (atom.attrib[kString]):
#                             kIndices.append(int(atom.attrib[kString]))
#                    except:
#                        pass
#
#                # set axis type based on k-Indices
#
#                charge = float(atom.attrib['c0'])
#
#                conversion = 1.0
#                dipole = [ conversion*float(atom.attrib['d1']), conversion*float(atom.attrib['d2']), conversion*float(atom.attrib['d3'])]
#
#                for t in types[0]:
#                    if (t not in generator.typeMap):
#                        generator.typeMap[t] = []
#
#                    valueMap = dict()
#                    valueMap['classIndex'] = atom.attrib['type']
#                    valueMap['kIndices'] = kIndices
#                    valueMap['charge'] = charge
#                    valueMap['dipole'] = dipole
#                    valueMap['quadrupole'] = quadrupole
#                    valueMap['axisType'] = axisType
#                    generator.typeMap[t].append(valueMap)
#
#            else:
#                outputString = "AmoebaMultipoleGenerator: error getting type for multipole: %s" % (atom.attrib['class'])
#                raise ValueError(outputString)
#
#        # polarization parameters
#
#        for atom in element.findall('Atom'):
#        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
#            types = forceField._findAtomTypes(atom, 1)
#            if None not in types:
#
#                classIndex = atom.attrib['type']
#                polarizability = float(atom.attrib['polarizability'])
#                charge = float(atom.attrib['charge'])
#                damping_factor = float(atom.attrib['damping-factor'])
#
#                for t in types[0]:
#                    if (t not in generator.typeMap):
#                        outputString = "AmoebaMultipoleGenerator: polarize type not present: %s" % (atom.attrib['type'])
#                        raise ValueError(outputString)
#                    else:
#                        typeMapList = generator.typeMap[t]
#                        hit = 0
#                        for (ii, typeMap) in enumerate(typeMapList):
#
#                            if (typeMap['classIndex'] == classIndex):
#                                typeMap['polarizability'] = polarizability
#                                typeMap['charge'] = charge
#                                typeMap['damping_factor'] = damping_factor
#                                typeMapList[ii] = typeMap
#                                hit = 1
#
#                        if (hit == 0):
#                            outputString = "MBPolElectrostaticsForceGenerator: error getting type for atom: class index=%s not in residue list?" % (atom.attrib['class'])
#                            raise ValueError(outputString)
#
#            else:
#                outputString = "MBPolElectrostaticsForceGenerator: error getting type for atom: %s" % (atom.attrib['class'])
#                raise ValueError(outputString)
#
#    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
#
#        methodMap = {NoCutoff:mm.AmoebaMultipoleForce.NoCutoff,
#                     PME:mm.AmoebaMultipoleForce.PME}
#
#        # get or create force depending on whether it has already been added to the system
#
#        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
#        existing = [f for f in existing if type(f) == mm.AmoebaMultipoleForce]
#        if len(existing) == 0:
#            force = mm.AmoebaMultipoleForce()
#            sys.addForce(force)
#            if (nonbondedMethod not in methodMap):
#                raise ValueError( "AmoebaMultipoleForce: input cutoff method not available." )
#            else:
#                force.setNonbondedMethod(methodMap[nonbondedMethod])
#            force.setCutoffDistance(nonbondedCutoff)
#
#            if ('ewaldErrorTolerance' in args):
#                force.setEwaldErrorTolerance(float(args['ewaldErrorTolerance']))
#
#            if ('polarization' in args):
#                polarizationType = args['polarization']
#                if (polarizationType.lower() == 'direct'):
#                    force.setPolarizationType(mm.AmoebaMultipoleForce.Direct)
#                else:
#                    force.setPolarizationType(mm.AmoebaMultipoleForce.Mutual)
#
#            if ('aEwald' in args):
#                force.setAEwald(float(args['aEwald']))
#
#            if ('pmeGridDimensions' in args):
#                force.setPmeGridDimensions(args['pmeGridDimensions'])
#
#            if ('mutualInducedMaxIterations' in args):
#                force.setMutualInducedMaxIterations(int(args['mutualInducedMaxIterations']))
#
#            if ('mutualInducedTargetEpsilon' in args):
#                force.setMutualInducedTargetEpsilon(float(args['mutualInducedTargetEpsilon']))
#
#        else:
#            force = existing[0]
#
#        # add particles to force
#        # throw error if particle type not available
#
#        # get 1-2, 1-3, 1-4, 1-5 bonded sets
#
#        # 1-2
#
#        bonded12ParticleSets = []
#        for i in range(len(data.atoms)):
#            bonded12ParticleSet = AmoebaVdwGenerator.getBondedParticleSet(i, data)
#            bonded12ParticleSet = set(sorted(bonded12ParticleSet))
#            bonded12ParticleSets.append(bonded12ParticleSet)
#
#        # 1-3
#
#        bonded13ParticleSets = []
#        for i in range(len(data.atoms)):
#            bonded13Set = set()
#            bonded12ParticleSet = bonded12ParticleSets[i]
#            for j in bonded12ParticleSet:
#                bonded13Set = bonded13Set.union(bonded12ParticleSets[j])
#
#            # remove 1-2 and self from set
#
#            bonded13Set = bonded13Set - bonded12ParticleSet
#            selfSet = set()
#            selfSet.add(i)
#            bonded13Set = bonded13Set - selfSet
#            bonded13Set = set(sorted(bonded13Set))
#            bonded13ParticleSets.append(bonded13Set)
#
#        # 1-4
#
#        bonded14ParticleSets = []
#        for i in range(len(data.atoms)):
#            bonded14Set = set()
#            bonded13ParticleSet = bonded13ParticleSets[i]
#            for j in bonded13ParticleSet:
#                bonded14Set = bonded14Set.union(bonded12ParticleSets[j])
#
#            # remove 1-3, 1-2 and self from set
#
#            bonded14Set = bonded14Set - bonded12ParticleSets[i]
#            bonded14Set = bonded14Set - bonded13ParticleSet
#            selfSet = set()
#            selfSet.add(i)
#            bonded14Set = bonded14Set - selfSet
#            bonded14Set = set(sorted(bonded14Set))
#            bonded14ParticleSets.append(bonded14Set)
#
#        # 1-5
#
#        bonded15ParticleSets = []
#        for i in range(len(data.atoms)):
#            bonded15Set = set()
#            bonded14ParticleSet = bonded14ParticleSets[i]
#            for j in bonded14ParticleSet:
#                bonded15Set = bonded15Set.union(bonded12ParticleSets[j])
#
#            # remove 1-4, 1-3, 1-2 and self from set
#
#            bonded15Set = bonded15Set - bonded12ParticleSets[i]
#            bonded15Set = bonded15Set - bonded13ParticleSets[i]
#            bonded15Set = bonded15Set - bonded14ParticleSet
#            selfSet = set()
#            selfSet.add(i)
#            bonded15Set = bonded15Set - selfSet
#            bonded15Set = set(sorted(bonded15Set))
#            bonded15ParticleSets.append(bonded15Set)
#
#        for (atomIndex, atom) in enumerate(data.atoms):
#            t = data.atomType[atom]
#            if t in self.typeMap:
#
#                multipoleList = self.typeMap[t]
#                hit = 0
#                savedMultipoleDict = 0
#
#                # assign multipole parameters via only 1-2 connected atoms
#
#                for multipoleDict in multipoleList:
#
#                    if (hit != 0):
#                        break
#
#                    kIndices = multipoleDict['kIndices']
#
#                    kz = kIndices[1]
#                    kx = kIndices[2]
#                    ky = kIndices[3]
#
#                    # assign multipole parameters
#                    #    (1) get bonded partners
#                    #    (2) match parameter types
#
#                    bondedAtomIndices = bonded12ParticleSets[atomIndex]
#                    zaxis = -1
#                    xaxis = -1
#                    yaxis = -1
#                    for bondedAtomZIndex in bondedAtomIndices:
#
#                       if (hit != 0):
#                           break
#
#                       bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
#                       bondedAtomZ = data.atoms[bondedAtomZIndex]
#                       if (bondedAtomZType == kz):
#                          for bondedAtomXIndex in bondedAtomIndices:
#                              if (bondedAtomXIndex == bondedAtomZIndex or hit != 0):
#                                  continue
#                              bondedAtomXType = int(data.atomType[data.atoms[bondedAtomXIndex]])
#                              if (bondedAtomXType == kx):
#                                  if (ky == 0):
#                                      zaxis = bondedAtomZIndex
#                                      xaxis = bondedAtomXIndex
#                                      if( bondedAtomXType == bondedAtomZType and xaxis < zaxis ):
#                                          swapI = zaxis
#                                          zaxis = xaxis
#                                          xaxis = swapI
#                                      else:
#                                          for bondedAtomXIndex in bondedAtomIndices:
#                                              bondedAtomX1Type = int(data.atomType[data.atoms[bondedAtomXIndex]])
#                                              if( bondedAtomX1Type == kx and bondedAtomXIndex != bondedAtomZIndex and bondedAtomXIndex < xaxis ):
#                                                  xaxis = bondedAtomXIndex
#
#                                      savedMultipoleDict = multipoleDict
#                                      hit = 1
#                                  else:
#                                      for bondedAtomYIndex in bondedAtomIndices:
#                                          if (bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit != 0):
#                                              continue
#                                          bondedAtomYType = int(data.atomType[data.atoms[bondedAtomYIndex]])
#                                          if (bondedAtomYType == ky):
#                                              zaxis = bondedAtomZIndex
#                                              xaxis = bondedAtomXIndex
#                                              yaxis = bondedAtomYIndex
#                                              savedMultipoleDict = multipoleDict
#                                              hit = 2
#
#                # assign multipole parameters via 1-2 and 1-3 connected atoms
#
#                for multipoleDict in multipoleList:
#
#                    if (hit != 0):
#                        break
#
#                    kIndices = multipoleDict['kIndices']
#
#                    kz = kIndices[1]
#                    kx = kIndices[2]
#                    ky = kIndices[3]
#
#                    # assign multipole parameters
#                    #    (1) get bonded partners
#                    #    (2) match parameter types
#
#                    bondedAtom12Indices = bonded12ParticleSets[atomIndex]
#                    bondedAtom13Indices = bonded13ParticleSets[atomIndex]
#
#                    zaxis = -1
#                    xaxis = -1
#                    yaxis = -1
#
#                    for bondedAtomZIndex in bondedAtom12Indices:
#
#                       if (hit != 0):
#                           break
#
#                       bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
#                       bondedAtomZ = data.atoms[bondedAtomZIndex]
#
#                       if (bondedAtomZType == kz):
#                          for bondedAtomXIndex in bondedAtom13Indices:
#
#                              if (bondedAtomXIndex == bondedAtomZIndex or hit != 0):
#                                  continue
#                              bondedAtomXType = int(data.atomType[data.atoms[bondedAtomXIndex]])
#                              if (bondedAtomXType == kx and bondedAtomZIndex in bonded12ParticleSets[bondedAtomXIndex]):
#                                  if (ky == 0):
#                                      zaxis = bondedAtomZIndex
#                                      xaxis = bondedAtomXIndex
#
#                                      # select xaxis w/ smallest index
#
#                                      for bondedAtomXIndex in bondedAtom13Indices:
#                                          bondedAtomX1Type = int(data.atomType[data.atoms[bondedAtomXIndex]])
#                                          if( bondedAtomX1Type == kx and bondedAtomXIndex != bondedAtomZIndex and bondedAtomZIndex in bonded12ParticleSets[bondedAtomXIndex] and bondedAtomXIndex < xaxis ):
#                                              xaxis = bondedAtomXIndex
#
#                                      savedMultipoleDict = multipoleDict
#                                      hit = 3
#                                  else:
#                                      for bondedAtomYIndex in bondedAtom13Indices:
#                                          if (bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit != 0):
#                                              continue
#                                          bondedAtomYType = int(data.atomType[data.atoms[bondedAtomYIndex]])
#                                          if (bondedAtomYType == ky and bondedAtomZIndex in bonded12ParticleSets[bondedAtomYIndex]):
#                                              zaxis = bondedAtomZIndex
#                                              xaxis = bondedAtomXIndex
#                                              yaxis = bondedAtomYIndex
#                                              savedMultipoleDict = multipoleDict
#                                              hit = 4
#
#                # assign multipole parameters via only a z-defining atom
#
#                for multipoleDict in multipoleList:
#
#                    if (hit != 0):
#                        break
#
#                    kIndices = multipoleDict['kIndices']
#
#                    kz = kIndices[1]
#                    kx = kIndices[2]
#
#                    zaxis = -1
#                    xaxis = -1
#                    yaxis = -1
#
#                    for bondedAtomZIndex in bondedAtom12Indices:
#
#                        if (hit != 0):
#                            break
#
#                        bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
#                        bondedAtomZ = data.atoms[bondedAtomZIndex]
#
#                        if (kx == 0 and kz == bondedAtomZType):
#                            kz = bondedAtomZIndex
#                            savedMultipoleDict = multipoleDict
#                            hit = 5
#
#                # assign multipole parameters via no connected atoms
#
#                for multipoleDict in multipoleList:
#
#                    if (hit != 0):
#                        break
#
#                    kIndices = multipoleDict['kIndices']
#
#                    kz = kIndices[1]
#
#                    zaxis = -1
#                    xaxis = -1
#                    yaxis = -1
#
#                    if (kz == 0):
#                        savedMultipoleDict = multipoleDict
#                        hit = 6
#
#                # add particle if there was a hit
#
#                if (hit != 0):
#
#                    atom.multipoleDict = savedMultipoleDict
#                    atom.polarizationGroups = dict()
#                    newIndex = force.addMultipole(savedMultipoleDict['charge'], savedMultipoleDict['dipole'], savedMultipoleDict['quadrupole'], savedMultipoleDict['axisType'],
#                                                                 zaxis, xaxis, yaxis, savedMultipoleDict['thole'], savedMultipoleDict['pdamp'], savedMultipoleDict['polarizability'])
#                else:
#                    raise ValueError("Atom %s of %s %d was not assigned." %(atom.name, atom.residue.name, atom.residue.index))
#            else:
#                raise ValueError('No multipole type for atom %s %s %d' % (atom.name, atom.residue.name, atom.residue.index))
#
#        # set polar groups
#
#        self.setPolarGroups(data, bonded12ParticleSets, force)
#
#parsers["MBPolElectrostaticsForce"] = MBPolElectrostaticsForceGenerator.parseElement

## @private
class MBPolOneBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolOneBodyForceGenerator()
        forceField._forces.append(generator)

        # <MBPolOneBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolOneBodyForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolOneBodyForce_template in element.findall('MBPolOneBodyForce'):
            types = forceField._findAtomTypes(MBPolOneBodyForce_template, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolOneBodyForce_template.attrib['class1'],
                                    MBPolOneBodyForce_template.attrib['class2'],
                                    MBPolOneBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolOneBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolOneBodyForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolOneBodyForce and match their name
                print([atom2.index, atom1.index, atom3.index])
                force.addOneBody(atom2.index, atom1.index, atom3.index)

forcefield.parsers["MBPolOneBodyForce"] = MBPolOneBodyForceGenerator.parseElement

## @private
class MBPolTwoBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolTwoBodyForceGenerator()
        forceField._forces.append(generator)

        # <MBPolTwoBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolTwoBodyForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolTwoBodyForce_template in element.findall('MBPolTwoBodyForce'):
            types = forceField._findAtomTypes(MBPolTwoBodyForce_template, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolTwoBodyForce_template.attrib['class1'],
                                    MBPolTwoBodyForce_template.attrib['class2'],
                                    MBPolTwoBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolTwoBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolTwoBodyForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolTwoBodyForce and match their name
                print([atom2.index, atom1.index, atom3.index])
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolTwoBodyForce"] = MBPolTwoBodyForceGenerator.parseElement

## @private
class MBPolThreeBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolThreeBodyForceGenerator()
        forceField._forces.append(generator)

        # <MBPolThreeBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolThreeBodyForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolThreeBodyForce_template in element.findall('MBPolThreeBodyForce'):
            types = forceField._findAtomTypes(MBPolThreeBodyForce_template, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolThreeBodyForce_template.attrib['class1'],
                                    MBPolThreeBodyForce_template.attrib['class2'],
                                    MBPolThreeBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolThreeBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolThreeBodyForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolThreeBodyForce and match their name
                print([atom2.index, atom1.index, atom3.index])
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolThreeBodyForce"] = MBPolThreeBodyForceGenerator.parseElement

## @private
class MBPolDispersionForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolDispersionForceGenerator()
        forceField._forces.append(generator)

        # <MBPolDispersionForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolDispersionForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolDispersionForce_template in element.findall('MBPolDispersionForce'):
            types = forceField._findAtomTypes(MBPolDispersionForce_template, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolDispersionForce_template.attrib['class1'],
                                    MBPolDispersionForce_template.attrib['class2'],
                                    MBPolDispersionForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolDispersionForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolDispersionForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolDispersionForce and match their name
                print([atom2.index, atom1.index, atom3.index])
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolDispersionForce"] = MBPolDispersionForceGenerator.parseElement
