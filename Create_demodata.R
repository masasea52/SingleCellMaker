# 必要なパッケージをインストール (もしインストールされていなければ)
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("utils", quietly = TRUE)) install.packages("utils") # zip関数用
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils") # gzip圧縮用
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr") # データ操作用
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble") # rownames_to_column用

# 1. パラメータ設定
n_cells <- 500 # 細胞の数をご要望通り500に設定
output_dir_name <- "demo_10x_data_immune_cells" # 出力フォルダ名
zip_file_name <- "demo_10x_data.zip" # ZIPファイル名 (アプリで読み込むため同じ名前に)

# 2. 実在する遺伝子名をロード
# 必ず、添付された 'genes.tsv' ファイルをRの現在の作業ディレクトリに置いてください
genes_df <- read.delim("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
  dplyr::rename(ID = V1, Symbol = V2)

# デモデータに使用する遺伝子の選択
n_total_genes_to_use <- 2000 # 使用する総遺伝子数
mt_genes_from_file <- genes_df$Symbol[grepl("^MT-", genes_df$Symbol, ignore.case = TRUE)]
non_mt_genes_from_file <- genes_df$Symbol[!grepl("^MT-", genes_df$Symbol, ignore.case = TRUE)]

# 細胞種特異的マーカー遺伝子のリスト
cell_type_markers <- list(
  "Neutrophil" = c("ELANE", "MPO", "FCGR3B", "CSF3R", "CEACAM8"),
  "Macrophage" = c("CD68", "CSF1R", "MRC1", "CD163", "ITGAM"), # ITGAM=CD11b
  "Lymphocyte" = c("CD3D", "CD3G", "CD79A", "MS4A1", "NKG7", "CD8A") # T, B, NK細胞マーカーを混ぜる
)
# all_marker_symbols は、今回のスクリプトでは直接は使われないが、参考として残す
# all_marker_symbols <- unique(unname(unlist(cell_type_markers)))

# 選択する遺伝子シンボル: マーカー遺伝子 + ランダムな非マーカー遺伝子 + MT遺伝子
selected_marker_symbols <- unique(unname(unlist(lapply(cell_type_markers, function(markers) {
  markers[markers %in% genes_df$Symbol] # genes.tsvに存在するマーカーのみ選択
}))))

# 残りの遺伝子枠を非マーカー遺伝子で埋める
num_random_genes_to_add <- n_total_genes_to_use - length(selected_marker_symbols) - min(5, length(mt_genes_from_file))
random_non_marker_genes <- setdiff(non_mt_genes_from_file, selected_marker_symbols)
if (length(random_non_marker_genes) < num_random_genes_to_add) {
  selected_random_symbols <- random_non_marker_genes
  num_mt_to_select <- n_total_genes_to_use - length(selected_marker_symbols) - length(selected_random_symbols)
  num_mt_to_select <- min(num_mt_to_select, length(mt_genes_from_file))
} else {
  selected_random_symbols <- sample(random_non_marker_genes, num_random_genes_to_add, replace = FALSE)
  num_mt_to_select <- min(5, length(mt_genes_from_file))
}

selected_gene_symbols_final <- unique(c(selected_marker_symbols, selected_random_symbols, sample(mt_genes_from_file, num_mt_to_select, replace = FALSE)))
# 選択された遺伝子シンボルに対応するIDを取得し、保存順を定義
selected_genes_df <- genes_df %>%
  dplyr::filter(Symbol %in% selected_gene_symbols_final) %>%
  dplyr::distinct(Symbol, .keep_all = TRUE) %>% # 同じ遺伝子シンボルがある場合に重複を避ける
  dplyr::arrange(Symbol) # アルファベット順に並べることで再現性を高める (または元のgenes.tsvの順序を保つ)

gene_symbols_final <- selected_genes_df$Symbol
gene_ids_final <- selected_genes_df$ID
n_genes_final <- length(gene_symbols_final)
message(paste0("Total genes selected: ", n_genes_final))
message(paste0("  - Marker genes (present in list): ", length(selected_marker_symbols)))
message(paste0("  - MT genes (present in list): ", sum(grepl("^MT-", gene_symbols_final))))


# 3. カウント行列の生成
set.seed(456) # カウント生成のシード
dummy_counts <- matrix(0, nrow = n_genes_final, ncol = n_cells)
rownames(dummy_counts) <- gene_symbols_final # 遺伝子シンボルを行名に設定
colnames(dummy_counts) <- paste0("cell_", 1:n_cells)

# 細胞を3つのクラスターに分割 (好中球, マクロファージ, リンパ球)
cell_types <- c("Neutrophil", "Macrophage", "Lymphocyte")
cell_clusters <- sample(0:2, n_cells, replace = TRUE, prob = c(0.3, 0.3, 0.4)) # 比率を調整
names(cell_clusters) <- colnames(dummy_counts)


# ベースの発現をランダムに割り当てる
for (i in 1:n_cells) {
  expressed_genes <- sample(1:n_genes_final, sample(50:200, 1)) # 各細胞で発現する遺伝子の数
  dummy_counts[expressed_genes, i] <- sample(1:20, length(expressed_genes), replace = TRUE)
}

# 細胞種特異的マーカー遺伝子の発現を強化
for (cell_type_idx in 0:2) {
  current_cell_type_name <- cell_types[cell_type_idx + 1]
  marker_genes_for_type <- cell_type_markers[[current_cell_type_name]]
  
  cells_in_current_cluster <- which(cell_clusters == cell_type_idx)
  
  if (length(cells_in_current_cluster) > 0 && length(marker_genes_for_type) > 0) {
    for (marker_symbol in marker_genes_for_type) {
      if (marker_symbol %in% rownames(dummy_counts)) { # 選択された遺伝子リストにマーカーが存在するか確認
        marker_row_idx <- which(rownames(dummy_counts) == marker_symbol)
        # 当該細胞種に属する細胞でマーカーを高発現させる
        dummy_counts[marker_row_idx, cells_in_current_cluster] <-
          pmax(dummy_counts[marker_row_idx, cells_in_current_cluster], 
               sample(20:50, length(cells_in_current_cluster), replace = TRUE))
        
        # 他の細胞種では発現を低く（またはゼロに）抑える
        cells_in_other_clusters <- setdiff(1:n_cells, cells_in_current_cluster)
        dummy_counts[marker_row_idx, cells_in_other_clusters] <-
          pmin(dummy_counts[marker_row_idx, cells_in_other_clusters], 
               sample(0:5, length(cells_in_other_clusters), replace = TRUE))
      }
    }
  }
}

# ミトコンドリア遺伝子の発現を意図的に高くする（例: リンパ球の一部で）
mt_gene_indices <- which(grepl("^MT-", rownames(dummy_counts), ignore.case = TRUE))
if (length(mt_gene_indices) > 0) {
  lymph_cells_for_high_mt <- sample(which(cell_clusters == 2), sum(cell_clusters == 2) / 2)
  if (length(lymph_cells_for_high_mt) > 0) {
    dummy_counts[mt_gene_indices, lymph_cells_for_high_mt] <-
      pmax(dummy_counts[mt_gene_indices, lymph_cells_for_high_mt], 
           sample(20:50, length(mt_gene_indices) * length(lymph_cells_for_high_mt), replace = TRUE))
  }
}

dummy_sparse_matrix <- as(dummy_counts, "dgCMatrix")

# 4. 出力ディレクトリの作成
output_dir <- file.path(getwd(), output_dir_name)
if (!dir.exists(output_dir)) dir.create(output_dir)

# 5. ファイルの保存
message("Saving files to: ", output_dir)

# a) matrix.mtx.gz (Matrix Market format, gzipped)
mtx_path <- file.path(output_dir, "matrix.mtx")
Matrix::writeMM(obj = dummy_sparse_matrix, file = mtx_path)
R.utils::gzip(filename = mtx_path, destname = paste0(mtx_path, ".gz"), remove = TRUE)

# b) features.tsv.gz
features_df_output <- data.frame(
  ID = gene_ids_final[match(rownames(dummy_sparse_matrix), gene_symbols_final)], # 保存された行列の行名順にIDを取得
  Name = rownames(dummy_sparse_matrix),
  Type = "Gene"
)
features_path <- file.path(output_dir, "features.tsv")
write.table(features_df_output,
            file = features_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
R.utils::gzip(filename = features_path, destname = paste0(features_path, ".gz"), remove = TRUE)

# c) barcodes.tsv.gz - **ここが簡潔化のポイント**
# cell_typeメタデータはSeuratオブジェクトに直接付与せず、barcodes.tsvには標準の1列のみを保存
barcodes_df_output <- data.frame(
  Barcode = colnames(dummy_sparse_matrix)
)
barcodes_path <- file.path(output_dir, "barcodes.tsv")
write.table(barcodes_df_output, # 標準の1列のみ保存
            file = barcodes_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
R.utils::gzip(filename = barcodes_path, destname = paste0(barcodes_path, ".gz"), remove = TRUE)

message("Successfully created 'matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz' in folder: ", output_dir)

# 6. 生成したフォルダをZIP圧縮
zip_path <- file.path(getwd(), zip_file_name)
files_to_zip <- list.files(output_dir, pattern = "\\.gz$", full.names = TRUE)
utils::zip(zipfile = zip_path, files = files_to_zip, flags = "-j")

message("Successfully created zip file: ", zip_path)
message("You can now copy this file to your Shiny app's 'www' folder.")

# 7. 生成した一時フォルダをクリーンアップ (任意)
# unlink(output_dir, recursive = TRUE)