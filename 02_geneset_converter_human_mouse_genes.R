# Install if needed
install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)

# Connect to Ensembl
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get human ortholog mappings
mapping <- getBM(
  attributes = c("mgi_symbol", "hsapiens_homolog_associated_gene_name"),
  mart = mouse
)

# Example mouse gene list
Shaul_DAM <- c("Apoe","Ctsd","Ctss","B2m","Gm10925","Ctsb","Tyrobp","Cst7",
               "Fth1","Gnas","Ctsl","Cd9","Grn","mmu-mir-703","Ctsz","Gm11410",
               "Trem2","Rplp1","Lyz2","AC151602.1","Mpeg1","Hexa")

# Convert to human gene names
Shaul_DAM_human <- mapping$hsapiens_homolog_associated_gene_name[
  match(Shaul_DAM, mapping$mgi_symbol)
]

# Combine into a data frame
converted <- data.frame(mouse_gene = Shaul_DAM, human_gene = Shaul_DAM_human)
print(converted)