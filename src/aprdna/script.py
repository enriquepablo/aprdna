

import argparse

from Bio.Seq import Seq
from reportlab.lib.units import cm
from Bio import SeqIO
from Bio.Graphics import BasicChromosome
from Bio.SeqFeature import SeqFeature, FeatureLocation


dinucleotides = ['AA', 'AT', 'AG', 'AC',
                 'TA', 'TT', 'TG', 'TC',
                 'GA', 'GT', 'GG', 'GC',
                 'CA', 'CT', 'CG', 'CC',
                 ]


class Search:

    def __init__(self, nucleotide, period, tolerance, min_length, record):
        self.nucleotide = Seq(nucleotide)
        self.period = period
        self.tolerance = tolerance
        self.min_length = min_length

        self.current_position = 0
        self.prev_position = 0
        self.current_search = []

        self.featured = 0

        self.record = record
        self.dna = self.record.seq
        self.length = len(self.dna)

        self.chr_diagram = BasicChromosome.Organism()
        self.chr_diagram.page_size = (29.7 * cm, 21 * cm)  # A4 landscape

        self.features = []

    def run(self):

        while self.current_position + self.period < self.length:
            found = self.find_start_position()
            if found:
                self.do_search()
            else:
                break

        per10000 = int((self.featured / self.length) * 10000)
        legend = f"{self.nucleotide}: {per10000:>5,d} %00 ({len(self.features):>5})"

        ftypes = set([f.type for f in self.record.features])
        for ftype in ftypes:
            features = []
            for f in self.record.features:
                if f.type == ftype:
                    f.qualifiers['color'] = [2]
                    features.append(f)
            self.add_chromosome(features, f"{ftype}")

        self.chr_diagram.draw(str(self.nucleotide) + '-' + self.record.name + '-' + str(self.min_length) + '-' + str(self.tolerance) + ".pdf", self.record.description)

        print(legend)

    def add_chromosome(self, other_feature, name):

        chromosome = BasicChromosome.Chromosome(name)
        chromosome.scale_num = self.length + 4

        start = BasicChromosome.TelomereSegment()
        start.scale = 2
        chromosome.add(start)

        features = self.features + other_feature

        body = BasicChromosome.AnnotatedChromosomeSegment(self.length, features)
        body.scale = self.length
        chromosome.add(body)

        # Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = 2
        chromosome.add(end)

        # This chromosome is done
        self.chr_diagram.add(chromosome)

    def find_start_position(self):
        found = False
        while not found:
            if self.current_position + self.period <= self.length:
                search_area = self.dna[self.current_position: self.current_position + self.period]

                try:
                    start_position = search_area.index(self.nucleotide)
                    found = True
                    self.current_position += start_position + 1
                    self.current_search.append(self.current_position)
                except ValueError:
                    self.current_position += self.period
            else:
                return False

        return True

    def do_search(self):
        finding = True
        while finding:
            finding, position = self.do_search_step()
            if finding:
                self.current_search.append(position)
                self.current_position = position

        if len(self.current_search) >= self.min_length:
            location = FeatureLocation(self.current_search[0], self.current_search[-1])
            feature = SeqFeature(location, type=self.nucleotide + " Island")
            self.features.append(feature)
            self.featured += self.current_search[-1] - self.current_search[0]

        self.current_search = []

    def do_search_step(self):
        next_start = int(self.current_position + self.period - (self.tolerance / 2))
        next_end = int(next_start + self.tolerance)

        search_area = self.dna[next_start: next_end]

        try:
            next_pos = search_area.index(self.nucleotide)
            return True, next_pos + self.current_position + 1

        except ValueError:
            return False, -1


def main():
    parser = argparse.ArgumentParser(description='Search for repetitions of nucleotides.')
    parser.add_argument('-p', '--period', help='Periodicity of the repetition', type=int)
    parser.add_argument('-t', '--tolerance', help='Tolerance window for the repetition', type=int)
    parser.add_argument('-m', '--length', help='Minimal repetition length', type=int)
    parser.add_argument('-r', '--throttle', help='Stop after so many fragments', type=int, default=-1)
    parser.add_argument('path', help='Path to Sequence file', type=str)

    args = parser.parse_args()

    print(args.path.split('/')[-1])

    for i, record in enumerate(SeqIO.parse(args.path, "genbank")):
        print('\n')
        print(f"{record.name} - {len(record.seq):,d} bp")
        if args.throttle != -1 and i >= args.throttle:
            break

        for nucleotide in dinucleotides:
            search = Search(nucleotide, args.period, args.tolerance, args.length, record)
            search.run()


if __name__ == "__main__":
    main()
