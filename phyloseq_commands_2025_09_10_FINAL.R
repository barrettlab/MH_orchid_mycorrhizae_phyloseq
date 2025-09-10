# Phyloseq

# Files needed:
#----1) otu_table_sintax.biom 			# OTU table in biom format
#----2) repseqs.fasta					# All unique OTU ITS2 sequences
#----3) metadata_w_reads.txt					# Tab-delimited metadata text file (sample, genusspecies, species, var, state, county, lat, lon)
#----4) repseqs_tree.nwk				# ***Newick tree file for repseqs. HOWEVER, I would wait to filter first, then export fasta, align, and build tree to re-import
										# Otherwise, you'll be aligning ~20,000 OTUs instead of a few thousand!

# Just a note, this dataset was generated from several species, so #### STEP 6 #### below is critical if you want to focus on one species.

library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(decontam)
library(FUNGuildR)
library(edgeR)

# Package citation info
# citation(phyloseq)
# citation(dplyr)
# citation(tidyr)
# citation(stringr)
# citation(ggplot2)
# citation(decontam)
# citation(FUNGuildR)
# citation(edgeR)

ps = import_biom("otu_table_sintax.biom", refseqfilename="repseqs.fasta")
meta_data <- import_qiime_sample_data("metadata_w_reads.txt")

## The ranks are as follows:

# Rank1 = Kingdom/Division
# Rank2 = Phylum
# Rank3 = Class
# Rank4 = Order
# Rank5 = Family
# Rank6 = Genus
# Rank7 = Species

## Merge the data and check the overview

ps2 <- merge_phyloseq(ps, meta_data)
sample_data(ps2)

# Function to use edgeR for TMM normalization: "trimmed mean of M values (TMM) that corrects the effective library size of the count tables"
tmm_normalize_phyloseq <- function(ps) {
  C <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) C <- t(C)
  Y <- edgeR::DGEList(counts = C)
  Y <- edgeR::calcNormFactors(Y, method = "TMM")
  eff <- Y$samples$lib.size * Y$samples$norm.factors
  list(
    eff_lib_size = eff,
    ps_tmm       = { tmp <- ps; otu_table(tmp) <- otu_table(t(t(C)/eff)*mean(eff), TRUE); tmp },
    ps_tmm_rel   = { tmp <- ps; otu_table(tmp) <- otu_table(t(t(C)/eff), TRUE); tmp },
    ps_tmm_cpm   = { tmp <- ps; otu_table(tmp) <- otu_table(edgeR::cpm(Y, TRUE), TRUE); tmp },
    ps_tmm_logcpm= { tmp <- ps; otu_table(tmp) <- otu_table(edgeR::cpm(Y, TRUE, log=TRUE, prior.count=1), TRUE); tmp }
  )
}

# usage:
res <- tmm_normalize_phyloseq(ps2)
ps_tmm      <- res$ps_tmm
ps_tmm_rel  <- res$ps_tmm_rel
eff.lib.size<- res$eff_lib_size

# 1) Keep only Fungi (handles Rank1..Rank7 with k__/p__/… prefixes)

tt_mat <- as.matrix(tax_table(ps_tmm))

find_rank_col <- function(tt, pattern) {
  hits <- sapply(colnames(tt), function(cl) {
    sum(grepl(pattern, tt[, cl], ignore.case = TRUE), na.rm = TRUE)
  })
  if (all(hits == 0)) NA_character_ else names(which.max(hits))
}

kcol <- find_rank_col(tt_mat, '^"?k__'); stopifnot(!is.na(kcol))
kingdom_raw <- tt_mat[, kcol, drop = TRUE]

fungal_taxa <- taxa_names(ps_tmm)[!is.na(kingdom_raw) &
                              grepl('^"?k__Fungi', kingdom_raw, ignore.case = TRUE)]
ps_fungi <- prune_taxa(fungal_taxa, ps_tmm)

cat("Start taxa:", ntaxa(ps_tmm), " -> Fungi-only:", ntaxa(ps_fungi), "\n")

# 2) Remove contaminants (decontam)

# BiocManager::install("decontam")

# If you have negatives or DNA conc, run decontam here on the appropriate object.
# Otherwise, carry forward fungi-only:
# ps_nocontam <- ps_fungi

# View library size distribution by species
df <- as.data.frame(sample_data(ps_fungi)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_fungi)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=species)) + geom_point()

# Find contaminants, using 'reads' as a proxy of concentration

contamdf.freq <- isContaminant(ps_fungi, method="frequency", conc="reads")
head(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

# Filter out contaminants detected by decontam
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps_fungi)
ps.noncontam



# 3) Light prevalence & abundance filtering

otu <- otu_table(ps.noncontam); if (!taxa_are_rows(ps.noncontam)) otu <- t(otu)
prev <- rowSums(otu > 0)
tot  <- rowSums(otu)

prev_min <- max(2, ceiling(0.02 * nsamples(ps.noncontam)))  # ≥2% of samples
tot_min  <- 10                                             # ≥10 reads total

keep <- prev >= prev_min & tot >= tot_min
ps_filt <- prune_taxa(keep, ps.noncontam)

cat("After filters:", ntaxa(ps_filt), "taxa\n")


# 4) Standardize taxonomy rank names (Rank1..Rank7 → Kingdom…Species)

tt <- as.matrix(tax_table(ps_filt))

find_rank_col <- function(tt, pattern) {
  hits <- sapply(colnames(tt), function(cl) sum(grepl(pattern, tt[, cl], ignore.case = TRUE), na.rm = TRUE))
  if (all(hits == 0)) NA_character_ else names(which.max(hits))
}
strip_rank_prefix <- function(x) {
  x <- as.character(x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^[a-z]__", "", x, ignore.case = TRUE)
  x[x %in% c("", "Unassigned", "unassigned")] <- NA
  x
}

rank_patterns <- c(Kingdom='^"?k__', Phylum='^"?p__', Class='^"?c__',
                   Order='^"?o__',  Family='^"?f__', Genus='^"?g__', Species='^"?s__')

new_names <- colnames(tt)
for (nm in names(rank_patterns)) {
  col <- find_rank_col(tt, rank_patterns[[nm]])
  if (!is.na(col)) new_names[new_names == col] <- nm
}
colnames(tt) <- new_names
for (cn in intersect(names(rank_patterns), colnames(tt))) tt[, cn] <- strip_rank_prefix(tt[, cn])

tax_table(ps_filt) <- tax_table(tt)

# Prefer Family, then Order, then Genus, etc., for plotting later
tax_cols_present <- colnames(tax_table(ps_filt))
priority <- c("Family","Order","Genus","Class","Phylum","Kingdom",
              "Species","Rank7","Rank6","Rank5","Rank4","Rank3","Rank2","Rank1")
fill_rank <- priority[priority %in% tax_cols_present][1]
fill_rank


# 5) Keep OMF (orchid mycorrhiza) + ECM (ectomycorrhizal) taxa

# Pull taxonomy as a data.frame for easy tagging
tx <- as.data.frame(tax_table(ps_filt), stringsAsFactors = FALSE)
if (!"Family" %in% names(tx)) tx$Family <- NA_character_
if (!"Order"  %in% names(tx)) tx$Order  <- NA_character_

# --- Orchid mycorrhiza (OMF) ---
omf_families <- c("Tulasnellaceae", "Ceratobasidiaceae", "Serendipitaceae", "Sebacinaceae")
omf_orders   <- c("Sebacinales", "Cantharellales")  # safety net if Family is missing

# --- Ectomycorrhiza (ECM) ---
ecm_families <- c(
  "Amanitaceae","Boletaceae","Suillaceae","Paxillaceae","Gomphidiaceae",
  "Rhizopogonaceae","Sclerodermataceae","Hymenogastraceae","Cortinariaceae",
  "Inocybaceae","Tricholomataceae","Hydnangiaceae","Cantharellaceae",
  "Clavulinaceae","Russulaceae","Thelephoraceae","Bankeraceae","Albatrellaceae",
  "Hydnaceae","Tuberaceae"
)

# --- Ericoid mycorrhiza (ERM) ---
# Most ERM taxa occur within Helotiales (Leotiomycetes); Oidiodendron is in Herpotrichiellaceae (Chaetothyriales).
erm_families <- c("Hyaloscyphaceae","Helotiaceae","Dermateaceae","Herpotrichiellaceae")
erm_orders   <- c("Helotiales")  # safety net (broad—may include non-ERM lineages)

# --- Arbuscular mycorrhiza (AM; Glomeromycota) ---
am_families <- c(
  "Glomeraceae","Acaulosporaceae","Gigasporaceae","Diversisporaceae",
  "Claroideoglomeraceae","Archaeosporaceae","Paraglomeraceae","Ambisporaceae",
  "Pacisporaceae","Geosiphonaceae"
)
am_orders <- c("Glomerales","Diversisporales","Paraglomerales","Archaeosporales","Ambisporales")  # safety net

# Union of family lists and order safety nets
myco_family_set <- unique(c(omf_families, ecm_families, erm_families, am_families))
myco_order_set  <- unique(c(omf_orders, erm_orders, am_orders))

# Normalize case for robust matching
canon <- function(x) stringr::str_to_title(trimws(as.character(x)))
fam  <- canon(tx$Family)
ord  <- canon(tx$Order)

is_myco <- (!is.na(fam) & fam %in% myco_family_set) |
           (!is.na(ord) & ord %in% myco_order_set)

# (Optional) store a label column; useful if you want to plot by guild later
mycotype <- dplyr::case_when(
  !is.na(fam) & fam %in% omf_families ~ "OMF",
  !is.na(fam) & fam %in% ecm_families ~ "ECM",
  !is.na(fam) & fam %in% erm_families ~ "ERM",
  !is.na(fam) & fam %in% am_families  ~ "AM",
  !is.na(ord) & ord %in% omf_orders   ~ "OMF (ord)",
  !is.na(ord) & ord %in% erm_orders   ~ "ERM (ord)",
  !is.na(ord) & ord %in% am_orders    ~ "AM (ord)",
  TRUE                                ~ "Other/unknown"
)

# Write Guild back (optional)
tax_table(ps_filt) <- tax_table(cbind(as.matrix(tax_table(ps_filt)), MycoType = mycotype))

# Keep only mycorrhizal taxa (family-based, with order safety net)
ps_myco <- subset_taxa(ps_filt, is_myco)
cat("After mycorrhizal keep-step (OMF+ECM+ERM+AM):", ntaxa(ps_myco), "taxa\n")

########################################################################################################################
########################################################################################################################
########################################################################################################################

# 6) (Optional) Subset to a focal ORCHID TAXON (e.g., genusspecies == "C-striata" # This only keeps Cephalanthera austiniae)

##### THIS IS WHERE YOU CHOOSE THE FOCAL TAXON #####

#### List the choices with:
unique(ps_myco@sam_data$genusspecies)

# [1] "C-aust"     "C-bent"     "C-involuta" "C-odont"    "C-striata"  "C-trif"    
# [7] "C-wist" 

# Comment out to use all hosts
ps_target <- subset_samples(ps_myco, genusspecies == "C-striata")   #### <<<< Critical line ####
ps_target <- prune_samples(sample_sums(ps_target) > 0, ps_target)
ps_target <- prune_taxa(taxa_sums(ps_target) > 0, ps_target)

# If you skipped the subset above, just do:
# ps_target <- ps_myco


########################################################################################################################
########################################################################################################################
########################################################################################################################



# 7) FUNGuildR 
Taxa_processed <- as.data.frame(tax_table(ps_target))
original_rownames <- rownames(Taxa_processed)
Taxa_processed <- Taxa_processed %>%
mutate(Taxonomy = paste0(Kingdom, ";", Phylum, ";", Class, ";", Order, ";",
Family, ";", Genus, ";", Species))
fung <- get_funguild_db()
Funguildresults <- FUNGuildR::funguild_assign(Taxa_processed, db = fung)

# Merge Funguild results back with phyloseq object:

# Ensure rownames match OTUs
rownames(Funguildresults) <- original_rownames

# Add to tax_table
tax_tab <- as.data.frame(tax_table(ps_target), stringsAsFactors = FALSE)
tax_tab$Guild <- Funguildresults$trophicMode

# guild instead (ugly!)
# tax_tab$Guild <- Funguildresults$guild

# Rebuild tax_table with guilds
tax_table(ps_target) <- tax_table(as.matrix(tax_tab))


# 8) Relative abundance and county barplot

ps_rel <- transform_sample_counts(ps_target, function(x) if (sum(x)==0) x else x / sum(x))

# # OPTIONAL: Filter out OTUs that never exceed 1% relative abundance in any sample
# ps_rel_1pct <- prune_taxa(
#     taxa_sums(ps_rel) > 0 & apply(otu_table(ps_rel), 1, function(x) any(x > 0.01)),
#     ps_rel
# )
# 
# # If you haven't already, get the relative abundances:
# 
# ps_rel <- transform_sample_counts(ps_rel_1pct, function(x) if (sum(x)==0) x else x / sum(x))

# (Optional) give NA families a label so they plot cleanly
tt <- as.data.frame(tax_table(ps_rel))
if ("Family" %in% names(tt)) tt$Family[is.na(tt$Family)] <- "Unassigned"
tax_table(ps_rel) <- tax_table(as.matrix(tt))



# 9) Saving rds for later analyses (remove '#' to uncomment)

#### Save the object
# saveRDS(ps_rel, file = "funguild_Caust_ps_rel.rds")

#### Load the object
# ps_rel <- readRDS("funguild_Caust_ps_rel.rds")

###################


# 10) Make a nice, interactive plot with plotly

# Step 1 — Add taxonomy + metadata into the plotting data frame
library(phyloseq)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(plotly)

# Package citation info
# citation(phyloseq)
# citation(ggplot2)
# citation(ggsci)
# citation(dplyr)
# citation(tidyr)
# citation(plotly)

# Extract data used for plotting
plot_data <- psmelt(ps_rel)  # Long-format data frame for plotting

# Combine taxonomy columns into one string
# (adjust based on what levels you have in tax_table)
plot_data <- plot_data %>%
  mutate(
    FullTaxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; "),
    FullTaxonomy = gsub("; NA", "", FullTaxonomy)  # clean missing values
  )

plot_data <- plot_data %>%
  mutate(ShortTaxonomy = paste(Family, Genus, Species, sep = "; "))


# Check the extra columns
head(plot_data[, c("OTU", "FullTaxonomy", "MycoType", "Guild")])

# With ggsci color scheme (IGV) to provide more contrast
# Change fille ="Family" to "Order" etc.

# Step 2 — Build the ggplot with custom tooltip aesthetics

p_county <- ggplot(plot_data, aes(
    x = sample,
    y = Abundance,
    fill = MycoType,
    # Pass tooltip fields to ggplotly via 'text' aesthetic
    text = paste0(
      "<b>Sample:</b> ", sample,
      "<br><b>Abundance:</b> ", round(Abundance, 3),
      "<br><b>Full taxonomy:</b> ", FullTaxonomy,
      "<br><b>MycoType:</b> ", MycoType,
      "<br><b>Guild:</b> ", Guild
    )
  )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_igv() +
  facet_wrap(~county, scales = "free_x") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    y = "Relative abundance",
    x = "Sample",
    title = "Corallorhiza striata samples — composition by Family (faceted by county)" # Change title to whatever taxon
  )

# Step 3 — Make it interactive with ggplotly()

pp <- ggplotly(p_county, tooltip = "text")
pp


