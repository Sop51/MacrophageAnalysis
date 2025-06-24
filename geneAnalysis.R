library(biomaRt)
library(ggplot2)
library(EnhancedVolcano)

# read in the datasets
salmonella_res <- read.csv('/Users/sm2949/Desktop/salmonella_results.csv', row.names = 1)
neural_res <- read.csv('/Users/sm2949/Desktop/neural_results.csv', row.names = 1)
heart_res <- read.csv('/Users/sm2949/Desktop/heart_results.csv', row.names = 1)

ensembl_symbol <- data.frame(
  ensembl_id = c(
    # Initial pro-inflammatory response
    "ENSDARG00000026925", "ENSDARG00000100791", "ENSDARG00000058557", "ENSDARG00000098700",
    
    # ECM remodelling enzymes and MMPs
    "ENSDARG00000044074", "ENSDARG00000001452", "ENSDARG00000002235", "ENSDARG00000008388",
    "ENSDARG00000043079", "ENSDARG00000077290", "ENSDARG00000010556", "ENSDARG00000051962",
    "ENSDARG00000013072", "ENSDARG00000033544",
    
    # Resolution of inflammation, growth and proliferation
    "ENSDARG00000100383", "ENSDARG00000078042", "ENSDARG00000027087", "ENSDARG00000104292",
    "ENSDARG00000086585", "ENSDARG00000075121", "ENSDARG00000031246",
    
    # Re-innervation, growth and proliferation
    "ENSDARG00000038609", "ENSDARG00000059387", "ENSDARG00000040982",
    
    # Pro-inflammatory and cytokine receptors
    "ENSDARG00000033727", "ENSDARG00000054542", "ENSDARG00000059294", "ENSDARG00000099902",
    "ENSDARG00000104314", "ENSDARG00000104474",
    
    # Repair and regeneration
    "ENSDARG00000099979", "ENSDARG00000061769",
    
    # ECM remodelling
    "ENSDARG00000044010", "ENSDARG00000100307", "ENSDARG00000019949"
  ),
  gene_symbol = c(
    # Initial pro-inflammatory response
    "nos2a", "tril", "il11b", "il16",
    
    # ECM remodelling enzymes and MMPs
    "loxl2b", "adam8a", "mmp14a", "mmp14b",
    "mmp23b", "mmp25a", "mmp25b", "mmp15a",
    "mmp15b", "adamts15b",
    
    # Resolution of inflammation
    "il10ra", "il10rb", "tgfb2", "ctgfb",
    "nrg2b", "hbegfa", "hbegfb",
    
    # Re-innervation
    "mpz", "fgf7", "fgf14",
    
    # Pro-inflammatory receptors
    "il12ba", "il12bb", "marco", "il17rc",
    "nrg1", "il6r",
    
    # Repair and regeneration
    "tgfbr3", "bmp10",
    
    # ECM remodelling
    "loxl2a", "adamts5", "serpinh1b"
  )
)


# pull out only these genes in heart model
heart_mac_cluster <- heart_res[rownames(heart_res) %in% ensembl_symbol$ensembl_id, ]
# add gene symbols based on Ensembl ID
heart_mac_cluster$gene_symbol <- ensembl_symbol$gene_symbol[match(rownames(heart_mac_cluster), ensembl_symbol$ensembl_id)]

EnhancedVolcano(
  heart_mac_cluster,
  lab = heart_mac_cluster$gene_symbol,           
  x = 'log2FoldChange',                           
  y = 'padj',                                    
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3,
  labSize = 4,
  colAlpha = 0.8,
  title = 'Heart Macrophage Injury Model',
  subtitle = 'Differential expression',
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5
)

# pull out only these genes in neural model
neural_mac_cluster <- neural_res[rownames(neural_res) %in% ensembl_symbol$ensembl_id, ]
# add gene symbols based on Ensembl ID
neural_mac_cluster$gene_symbol <- ensembl_symbol$gene_symbol[match(rownames(neural_mac_cluster), ensembl_symbol$ensembl_id)]

EnhancedVolcano(
  neural_mac_cluster,
  lab = neural_mac_cluster$gene_symbol,           
  x = 'logFC..Group..Injury.vs.Sham.',                           
  y = 'FDR..Group..Injury.vs.Sham.',                                    
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3,
  labSize = 4,
  colAlpha = 0.8,
  title = 'Neural Macrophage Injury Model',
  subtitle = 'Differential expression',
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5
)

# pull out only these genes in salmonella model
salmonella_mac_cluster <- salmonella_res[rownames(salmonella_res) %in% ensembl_symbol$gene_symbol, ]

EnhancedVolcano(
  salmonella_mac_cluster,
  lab = rownames(salmonella_mac_cluster),           
  x = 'log2FoldChange',                           
  y = 'padj',                                 
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3,
  labSize = 4,
  colAlpha = 0.8,
  title = 'Salmonella Macrophage Injury Model',
  subtitle = 'Differential expression',
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5
)
