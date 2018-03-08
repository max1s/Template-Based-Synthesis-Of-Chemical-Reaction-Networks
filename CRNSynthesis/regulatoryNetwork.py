from CRNSynthesis.symbolicLNA import *

class Network:
    def __init__(self):
        self.edges = []
        self.species = set()

    def add_edge(self, source, target, interaction_type):
        self.species.add(source)
        self.species.add(target)

        self.edges.append((source, target, interaction_type))

    def group_interactions(self):
        activating_species = {}  # dict in which each value is a list of species that activate the corresponding key
        repressing_species = {}
        for edge in self.edges:
            source, target, interaction_type = edge

            if interaction_type == "activate":
                a = activating_species
            else:
                a = repressing_species

            if target not in a.keys():
                a[target] = []
            a[target].append(source)

        return activating_species, repressing_species

    def to_crn(self, include_mrna=False):

        activating_species, repressing_species = self.group_interactions()

        # add species
        species_objects = {}
        for species_name in self.species:
            species_objects[species_name] = Species(species_name)

            if include_mrna:
                species_objects[species_name + "_mrna"] = Species(species_name + "_mrna")

        # add reactions
        reactions = []
        constant_offset = 1
        for species_name in self.species:

            activator_term = "0"
            activator_constants = []
            if species_name in activating_species.keys():
                activators = activating_species[species_name]
                activator_constants = ["k_%s" % (i + constant_offset) for i in range(len(activators))]
                activator_term = " + ".join(["%s * %s^2" % x for x in zip(activator_constants, activators)])
                constant_offset += len(activators)

            repressor_term = "0"
            repressor_constants = []
            if species_name in repressing_species.keys():
                repressors = repressing_species[species_name]
                repressor_constants = ["k_%s" % (i + constant_offset) for i in range(len(repressors))]
                repressor_term = " + ".join(["%s * %s^2" % x for x in zip(repressor_constants, repressors)])
                constant_offset += len(repressors)

            numerator = "k_%s + %s" % (species_name, activator_term)
            denominator = "1 + k_%s + %s + %s" % (species_name, activator_term, repressor_term)
            rate = "R_%s * (%s) / (%s)" % (species_name, numerator, denominator)

            all_constants = ["R_" + species_name, "k_" + species_name]
            all_constants.extend(activator_constants)
            all_constants.extend(repressor_constants)
            constant_objects = [RateConstant(x, 0, 1) for x in all_constants]

            if include_mrna:
                mrna = species_objects[species_name + "_mrna"]
                protein = species_objects[species_name]

                transcription = ArbitraryRateReaction([], [(mrna, 1)], rate, constant_objects)
                rna_degredation = Reaction([(mrna, 1)], [], RateConstant('d_' + species_name + "_mrna", 0, 1))
                translation = Reaction([(mrna, 1)], [(protein, 1), (mrna, 1)], RateConstant('tl_' + species_name, 0, 1))
                protein_degredation = Reaction([(protein, 1)], [], RateConstant('d_' + species_name, 0, 1))

                reactions.extend([transcription, rna_degredation, translation, protein_degredation])

            else:
                protein = species_objects[species_name]

                transcription_translation = ArbitraryRateReaction([], [(protein, 1)], rate, constant_objects)
                protein_degredation = Reaction([(protein, 1)], [], RateConstant('d_' + species_name, 0, 1))
                reactions.extend([transcription_translation, protein_degredation])

        return CRNSketch(reactions, [], [])

