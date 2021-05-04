import os
import argparse as argv

parser = argv.ArgumentParser(description="The input fastq folder that contains filtered fastq files for crispresso to work on")
parser.add_argument("Input_folder")
parser.add_argument("Ref_file")
args = parser.parse_args()
input_files = os.listdir(args.Input_folder)
amplicon_file = args.Ref_file
if args.Input_folder.endswith("/") == True:
	pass
else:
	args.Input_folder += "/"

out_dir = "/".join(args.Input_folder.split("/")[:-2])
# print(out_dir)
# print(args.Input_folder)
# print(input_files)

target_sites = []

with open(amplicon_file, "r") as f:
	lines = f.readlines()
	for line in lines:
		if line.startswith("TS"):
			target_sites.append(line.split(":")[1].rstrip("\n"))
		elif line.startswith("Name"):
			amplicon_name = line.split(":")[1].rstrip("\n")
		elif line.startswith("Sequence"):
			amplicon_seq = line.split(":")[1].rstrip("\n")


# print(amplicon_name)
# print(amplicon_seq)
# print(target_sites)

qwc = "_".join(target_sites)
print(qwc)

for fastq_file in input_files:
	stem = fastq_file.split("_")[0]
	command = "CRISPResso --fastq_r1 {} --amplicon_seq {} -o {} --fastq_output -n {} -qwc {} --ignore_substitutions".format(fastq_file, amplicon_seq, out_dir, stem, qwc)
	#print(command)
