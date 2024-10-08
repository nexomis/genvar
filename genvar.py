#!/usr/bin/env python

import argparse
import sys
import random
import numpy as np
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
  parser = argparse.ArgumentParser(description="Simulate viral evolution with mutations over passages")
  parser.add_argument("virus_fasta", help="Reference virus sequence in FASTA format")
  parser.add_argument("config_file", help="Configuration file (YAML format)")
  parser.add_argument("--num_fragments", type=int, default=100000, help="Number of fragments to simulate (default 1000000)")
  parser.add_argument("--R1_size", type=int, default=150, help="Read 1 size (default 150)")
  parser.add_argument("--R2_size", type=int, default=150, help="Read 2 size (default 150)")
  parser.add_argument("--fragment_size_mean", type=float, default=500, help="Mean fragment size (default 350)")
  parser.add_argument("--fragment_size_sd", type=float, default=20, help="Fragment size standard deviation (default 50)")
  parser.add_argument("--fragment_size_min", type=float, default=250, help="Fragment size standard deviation (default 250)")
  args = parser.parse_args()

  # Read virus fasta
  reference_seq_record = SeqIO.read(args.virus_fasta, "fasta")
  reference_sequence = str(reference_seq_record.seq)
  virus_length = len(reference_sequence)

  # Read configuration file
  with open(args.config_file, 'r') as f:
      config = yaml.safe_load(f)

  # Get passages from config
  passages = sorted(config.keys())

  # Initialize known mutations
  known_snps = []
  known_indels = []
  known_deletions = []

  for passage in passages:
      print(f"Processing passage: {passage}")
      passage_params = config[passage]
      # Generate new mutations
      snp_params = passage_params.get('snp', {})
      indel_params = passage_params.get('indel', {})
      deletion_params = passage_params.get('deletion', {})

      # Generate new mutations
      new_snps = generate_snp(snp_params, reference_sequence, known_snps, known_indels, known_deletions) if snp_params else []
      new_indels = generate_indel(indel_params, reference_sequence, known_snps, known_indels, known_deletions) if indel_params else []
      new_deletions = generate_deletions(deletion_params, reference_sequence, known_snps, known_indels, known_deletions) if deletion_params else []

      # Update existing mutations (mutation propagation)
      # For simplicity, we'll assume a Bernoulli trial with p=0.5 to determine if a mutation evolves
      mutation_stability_prob = 0.5

      for mutation_list, params in [(known_snps, snp_params), (known_indels, indel_params), (known_deletions, deletion_params)]:
          for mutation in mutation_list:
              if random.random() < mutation_stability_prob:
                  # Mutation remains stable
                  pass
              else:
                  # Mutation evolves
                  if params:
                      mutation.proportion = np.random.beta(params['beta_prop']['alpha'], params['beta_prop']['beta'])
                  else:
                      # If params not defined for this passage, mutation remains stable
                      pass

      # Add new mutations to known mutations
      known_snps.extend(new_snps)
      known_indels.extend(new_indels)
      known_deletions.extend(new_deletions)

      # Generate reads
      reads_R1, reads_R2 = generate_reads(reference_sequence, args.num_fragments, args.R1_size, args.R2_size, known_snps, known_indels, known_deletions, args.fragment_size_mean, args.fragment_size_sd, args.fragment_size_min)

      # Write mutations details to CSV
      with open(f"{passage}.snp.csv", 'w') as f_snp:
          f_snp.write("position,alt_base,proportion\n")
          for snp in known_snps:
              f_snp.write(f"{snp.position},{snp.alt_base},{snp.proportion}\n")

      with open(f"{passage}.indel.csv", 'w') as f_indel:
          f_indel.write("position,type,sequence,proportion\n")
          for indel in known_indels:
              f_indel.write(f"{indel.position},{indel.indel_type},{indel.sequence},{indel.proportion}\n")

      with open(f"{passage}.deletion.csv", 'w') as f_del:
          f_del.write("start_position,length,proportion\n")
          for deletion in known_deletions:
              f_del.write(f"{deletion.start_position},{deletion.length},{deletion.proportion}\n")

      # Write reads to FASTQ files
      R1_filename = f"{passage}_R1.fq.gz"
      R2_filename = f"{passage}_R2.fq.gz"

      with open(R1_filename, 'w') as f_R1:
          SeqIO.write(reads_R1, f_R1, "fastq")

      with open(R2_filename, 'w') as f_R2:
          SeqIO.write(reads_R2, f_R2, "fastq")

      print(f"Passage {passage} processing complete.")

class SNP:
  def __init__(self, position, alt_base, proportion):
      self.position = position  # position in reference, 0-based
      self.alt_base = alt_base
      self.proportion = proportion

  def affect_fragment(self, frag_start, frag_end):
      return self.position >= frag_start and self.position < frag_end

  def apply(self, frag_seq, frag_start):
      pos_in_frag = self.position - frag_start
      frag_seq = frag_seq[:pos_in_frag] + self.alt_base + frag_seq[pos_in_frag+1:]
      return frag_seq

  def overlap(self, pos):
      return self.position == pos

class Indel:
  def __init__(self, position, indel_type, sequence, proportion):
      self.position = position  # position in reference, 0-based
      self.indel_type = indel_type  # 'insertion' or 'deletion'
      self.sequence = sequence  # for insertion, the inserted seq; for deletion, sequence being deleted
      self.proportion = proportion
      self.length = len(sequence)

  def affect_fragment(self, frag_start, frag_end):
      if self.indel_type == 'insertion':
          return self.position >= frag_start and self.position <= frag_end  # insertion between bases
      else:
          del_end = self.position + self.length
          return (self.position >= frag_start and self.position < frag_end) or (del_end > frag_start and del_end <= frag_end)

  def apply(self, frag_seq, frag_start):
      pos_in_frag = self.position - frag_start
      if self.indel_type == 'insertion':
          frag_seq = frag_seq[:pos_in_frag] + self.sequence + frag_seq[pos_in_frag:]
      else:
          del_len = self.length
          frag_seq = frag_seq[:pos_in_frag] + frag_seq[pos_in_frag + del_len:]
      return frag_seq

  def overlap(self, pos):
      if self.indel_type == 'insertion':
          return self.position == pos
      else:
          return pos >= self.position and pos < self.position + self.length

class Deletion:
  def __init__(self, start_position, length, proportion):
      self.start_position = start_position  # position in reference, 0-based
      self.length = length
      self.proportion = proportion

  def affect_fragment(self, frag_start, frag_end):
      del_end = self.start_position + self.length
      return (self.start_position >= frag_start and self.start_position < frag_end) or (del_end > frag_start and del_end <= frag_end)

  def apply(self, frag_seq, frag_start):
      pos_in_frag = max(0, self.start_position - frag_start)
      del_start_in_frag = pos_in_frag
      del_end_in_frag = min(len(frag_seq), self.start_position + self.length - frag_start)
      frag_seq = frag_seq[:del_start_in_frag] + frag_seq[del_end_in_frag:]
      return frag_seq

  def overlap(self, pos):
      return pos >= self.start_position and pos < self.start_position + self.length

def generate_snp(snp_params, reference_sequence, known_snp=[], known_indel=[], known_deletions=[]):
  num_snps = int(np.random.gamma(shape=snp_params['gamma_pos']['shape'], scale=1.0/snp_params['gamma_pos']['rate']))
  positions = []
  snps = []
  # Generate positions ensuring they don't overlap with known mutations
  possible_positions = set(range(len(reference_sequence)))

  # Exclude positions of known mutations
  for snp in known_snp:
      possible_positions.discard(snp.position)
  for indel in known_indel:
      if indel.indel_type == 'insertion':
          possible_positions.discard(indel.position)
      else:
          for pos in range(indel.position, indel.position + indel.length):
              possible_positions.discard(pos)
  for deletion in known_deletions:
      for pos in range(deletion.start_position, deletion.start_position + deletion.length):
          possible_positions.discard(pos)

  positions = random.sample(tuple(possible_positions), min(num_snps, len(possible_positions)))

  for pos in positions:
      ref_base = reference_sequence[pos]
      alt_bases = [b for b in ['A', 'T', 'C', 'G'] if b != ref_base]
      alt_base = random.choice(alt_bases)
      proportion = np.random.beta(snp_params['beta_prop']['alpha'], snp_params['beta_prop']['beta'])
      snp = SNP(position=pos, alt_base=alt_base, proportion=proportion)
      snps.append(snp)
  return snps

def generate_indel(indel_params, reference_sequence, known_snp=[], known_indel=[], known_deletions=[]):
  num_indels = int(np.random.gamma(shape=indel_params['gamma_pos']['shape'], scale=1.0/indel_params['gamma_pos']['rate']))
  indels = []

  # Exclude positions of known mutations
  possible_positions = set(range(len(reference_sequence)))
  for snp in known_snp:
      possible_positions.discard(snp.position)
  for indel in known_indel:
      if indel.indel_type == 'insertion':
          possible_positions.discard(indel.position)
      else:
          for pos in range(indel.position, indel.position + indel.length):
              possible_positions.discard(pos)
  for deletion in known_deletions:
      for pos in range(deletion.start_position, deletion.start_position + deletion.length):
          possible_positions.discard(pos)

  positions = random.sample(tuple(possible_positions), min(num_indels, len(possible_positions)))

  for pos in positions:
      indel_type = 'insertion' if np.random.rand() < indel_params['bernouilli']['p'] else 'deletion'
      indel_size = random.choices([1,2,3], weights=indel_params['weights'])[0]
      proportion = np.random.beta(indel_params['beta_prop']['alpha'], indel_params['beta_prop']['beta'])
      if indel_type == 'insertion':
          # Create a random insertion sequence
          ins_seq = ''.join(random.choices(['A','T','C','G'], k=indel_size))
          indel = Indel(position=pos, indel_type='insertion', sequence=ins_seq, proportion=proportion)
      else:
          # For deletion, ensure the deletion does not go beyond the sequence
          if pos + indel_size > len(reference_sequence):
              indel_size = len(reference_sequence) - pos
          del_seq = reference_sequence[pos:pos+indel_size]
          indel = Indel(position=pos, indel_type='deletion', sequence=del_seq, proportion=proportion)
      indels.append(indel)
  return indels

def generate_deletions(deletion_params, reference_sequence, known_snp=[], known_indel=[], known_deletions=[]):
  num_deletions = int(np.random.gamma(shape=deletion_params['gamma_pos']['shape'], scale=1.0/deletion_params['gamma_pos']['rate']))
  deletions = []
  # Exclude positions of known mutations
  possible_positions = set(range(len(reference_sequence)))
  for snp in known_snp:
      possible_positions.discard(snp.position)
  for indel in known_indel:
      if indel.indel_type == 'insertion':
          possible_positions.discard(indel.position)
      else:
          for pos in range(indel.position, indel.position + indel.length):
              possible_positions.discard(pos)
  for deletion in known_deletions:
      for pos in range(deletion.start_position, deletion.start_position + deletion.length):
          possible_positions.discard(pos)

  positions = random.sample(tuple(possible_positions), min(num_deletions, len(possible_positions)))

  for pos in positions:
      del_size = int(np.random.gamma(shape=deletion_params['gamma_size']['shape'], scale=1.0/deletion_params['gamma_size']['rate']))
      if pos + del_size > len(reference_sequence):
          del_size = len(reference_sequence) - pos
      proportion = np.random.beta(deletion_params['beta_prop']['alpha'], deletion_params['beta_prop']['beta'])
      deletion = Deletion(start_position=pos, length=del_size, proportion=proportion)
      deletions.append(deletion)
  return deletions

def generate_frag_seq(reference_sequence, frag_start, frag_size, known_snps, known_indels, known_deletions):
  frag_seq = reference_sequence[frag_start:frag_start+frag_size]
  # Collect mutations that affect this fragment
  mutations = []
  # Sort mutations by position descending
  mutations.extend([snp for snp in known_snps if snp.affect_fragment(frag_start, frag_start+frag_size)])
  mutations.extend([indel for indel in known_indels if indel.affect_fragment(frag_start, frag_start+frag_size)])
  mutations.extend([deletion for deletion in known_deletions if deletion.affect_fragment(frag_start, frag_start+frag_size)])

  # Apply deletions and indels from highest position to lowest
  # For sorting, get position attribute (snp.position), (indel.position), (deletion.start_position)
  def get_mutation_position(mutation):
      if isinstance(mutation, SNP):
          return mutation.position
      elif isinstance(mutation, Indel):
          return mutation.position
      elif isinstance(mutation, Deletion):
          return mutation.start_position
      else:
          return 0

  mutations.sort(key=lambda x: get_mutation_position(x), reverse=True)
  for mutation in mutations:
      # Decide whether to apply mutation based on its proportion
      if np.random.rand() < mutation.proportion:
          frag_seq = mutation.apply(frag_seq, frag_start)
  return frag_seq

def generate_reads(reference_sequence, num_fragments, R1_size, R2_size, known_snps, known_indels, known_deletions, size_mean, size_sd, min_size):
  virus_length = len(reference_sequence)
  reads_R1 = []
  reads_R2 = []
  for i in range(num_fragments):
      # Generate fragment start position
      frag_size = 0
      while frag_size > min_size:
        frag_size = int(np.random.normal(size_mean, size_sd))
      if frag_size <= 0:
          frag_size = int(size_mean)
      if frag_size > virus_length:
          frag_size = virus_length
      frag_start = np.random.randint(0, virus_length - frag_size + 1)
      frag_seq = generate_frag_seq(reference_sequence, frag_start, frag_size, known_snps, known_indels, known_deletions)

      # Generate R1 and R2 sequences
      R1_seq = frag_seq[:R1_size]
      R2_seq = frag_seq[-R2_size:]  # Last R2_size bases
      # Reverse complement R2
      R2_seq = str(Seq(R2_seq).reverse_complement())

      # Generate quality scores (Phred score 28)
      R1_qual = [28]*len(R1_seq)
      R2_qual = [28]*len(R2_seq)

      # Create SeqRecord objects
      read_id = f"frag_{i}"
      read_R1 = SeqRecord(Seq(R1_seq), id=read_id + "/1", description="")
      read_R1.letter_annotations["phred_quality"] = R1_qual
      read_R2 = SeqRecord(Seq(R2_seq), id=read_id + "/2", description="")
      read_R2.letter_annotations["phred_quality"] = R2_qual

      reads_R1.append(read_R1)
      reads_R2.append(read_R2)
  return reads_R1, reads_R2

if __name__ == "__main__":
  main()
