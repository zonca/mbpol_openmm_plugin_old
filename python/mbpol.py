import mbpolplugin
from simtk.openmm.app import forcefield
from simtk import unit

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
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolDispersionForce"] = MBPolDispersionForceGenerator.parseElement
