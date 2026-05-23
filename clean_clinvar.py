import re

input_vcf = "clinvar.vcf"
output_file = "clinvar_simp.txt"

count_total = 0
count_omim = 0
count_pathogenic = 0
count_snv = 0

with open(input_vcf, "r", encoding="utf-8") as infile, open(output_file, "w", encoding="utf-8") as outfile:
    outfile.write("Name\tGene\tChromosome\tStart\tEnd\n")

    for line in infile:
        if line.startswith("#"):
            continue

        count_total += 1

        cols = line.strip().split("\t")

        chrom = cols[0]
        pos = cols[1]
        ref = cols[3]
        alt = cols[4]
        info = cols[7]

        # filtro 1: asociación con OMIM
        if "OMIM:" not in info:
            continue
        count_omim += 1

        # filtro 2: pathogenic solamente
        if "CLNSIG=Pathogenic" not in info:
            continue
        count_pathogenic += 1

        # filtro 3: solo SNVs
        if len(ref) != 1 or len(alt) != 1:
            continue
        count_snv += 1

        # extraer nombre HGVS
        hgvs_match = re.search(r'CLNHGVS=([^;]+)', info)
        name = hgvs_match.group(1) if hgvs_match else "NA"

        # extraer gen
        gene_match = re.search(r'GENEINFO=([^;]+)', info)
        gene = gene_match.group(1).split(":")[0] if gene_match else "NA"

        outfile.write(f"{name}\t{gene}\t{chrom}\t{pos}\t{pos}\n")

print("Filtrado terminado")
print("Total variantes:", count_total)
print("Con OMIM:", count_omim)
print("Pathogenic:", count_pathogenic)
print("SNVs finales:", count_snv)