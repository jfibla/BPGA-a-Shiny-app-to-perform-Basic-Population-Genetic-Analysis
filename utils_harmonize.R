# ─────────────────────────────────────────────
# utils_harmonize.R — Harmonització d'al·lels entre USER.bim i Reference.bim
# ─────────────────────────────────────────────

message("✅ utils_harmonize.R loaded")

suppressPackageStartupMessages({
  library(data.table)
})

# Complement d'al·lels
.comp_map <- c(A="T", T="A", C="G", G="C", "0"="0")

# Al·lels palindròmics (A/T o C/G, sin importar orden)
.is_palindromic <- function(a1, a2) {
  alleles <- paste(pmin(a1, a2), pmax(a1, a2))
  alleles %chin% c("A T","T A","C G","G C")
}

#' Harmoniza USER.bim contra Reference.bim y escribe out_bim_path (sin extensión)
#' El .bim de salida queda en el MISMO orden que USER.bim (filtrado),
#' para que coincida con el orden de PLINK al hacer --extract.
#'
#' @param user_bim_path Ruta a USER.bim
#' @param ref_bim_path  Ruta a Reference.bim
#' @param out_bim_path  Ruta base de salida (sin .bim)
#' @param drop_palindromic logical: descartar palindrómicos (por defecto TRUE)
#' @param dedup_strategy "first" o "error": cómo tratar SNP duplicados
#' @return lista invisible con n_input, n_kept, harmonized (data.table) y removed (data.table)
harmonize_bim_to_reference <- function(
    user_bim_path,
    ref_bim_path,
    out_bim_path,
    drop_palindromic = TRUE,
    dedup_strategy = c("first","error")[1]
) {
  # --- Carga segura ---
  if (!file.exists(user_bim_path)) stop("❌ USER.bim no encontrado: ", user_bim_path)
  if (!file.exists(ref_bim_path))  stop("❌ Reference.bim no encontrado: ", ref_bim_path)
  
  usr <- fread(user_bim_path, header = FALSE, na.strings = c("","NA"))
  if (ncol(usr) < 6) stop("❌ USER.bim con columnas insuficientes")
  setnames(usr, c("CHR","SNP","CM","BP","A1.s","A2.s"))
  
  ref <- fread(ref_bim_path, header = FALSE, na.strings = c("","NA"))
  if (ncol(ref) < 6) stop("❌ Reference.bim con columnas insuficientes")
  setnames(ref, c("CHR","SNP","CM","BP","A1.r","A2.r"))
  
  # Tipos básicos
  suppressWarnings({
    usr[, `:=`(CHR = as.character(CHR), BP = as.integer(BP), CM = as.numeric(CM))]
    ref[, `:=`(CHR = as.character(CHR), BP = as.integer(BP), CM = as.numeric(CM))]
  })
  
  # Normaliza alelos a mayúsculas
  for (j in c("A1.s","A2.s")) set(usr, j = j, value = toupper(usr[[j]]))
  for (j in c("A1.r","A2.r")) set(ref, j = j, value = toupper(ref[[j]]))
  
  # Validación de alelos (solo A,C,G,T,0)
  valid_alleles <- names(.comp_map)
  bad_usr <- which(!(usr$A1.s %in% valid_alleles & usr$A2.s %in% valid_alleles))
  if (length(bad_usr)) {
    warning("⚠️ Alelos no reconocidos en USER.bim (se descartarán en la fusión): filas ", paste(head(bad_usr,5), collapse=", "))
  }
  bad_ref <- which(!(ref$A1.r %in% valid_alleles & ref$A2.r %in% valid_alleles))
  if (length(bad_ref)) {
    warning("⚠️ Alelos no reconocidos en Reference.bim: filas ", paste(head(bad_ref,5), collapse=", "))
  }
  
  # Duplicados
  dup_usr <- usr[duplicated(SNP), SNP]
  dup_ref <- ref[duplicated(SNP), SNP]
  if (length(dup_usr) || length(dup_ref)) {
    msg <- sprintf("⚠️ Duplicados — USER: %d, REF: %d", length(dup_usr), length(dup_ref))
    if (identical(dedup_strategy, "error")) stop(msg) else warning(msg)
    if (dedup_strategy == "first") {
      usr <- usr[!duplicated(SNP)]
      ref <- ref[!duplicated(SNP)]
    }
  }
  
  # Índice para preservar orden de usuario
  usr[, idx := .I]
  
  setkey(usr, SNP)
  setkey(ref, SNP)
  
  # Join por SNP (inner)
  m <- ref[usr, nomatch = 0L]  # columnas de ref + usr alineadas por SNP
  
  # Renombra columnas para claridad post join
  # data.table mantiene columnas: de ref primero, luego de usr (sin duplicar 'SNP')
  # Asegúrate de tener A1.r, A2.r, A1.s, A2.s, idx
  if (!all(c("A1.r","A2.r","A1.s","A2.s","idx") %in% names(m)))
    stop("❌ Columnas esperadas no presentes tras el join.")
  
  # Sets de alelos (no orientados)
  m[, set_usr := paste(pmin(A1.s, A2.s), pmax(A1.s, A2.s))]
  m[, set_ref := paste(pmin(A1.r, A2.r), pmax(A1.r, A2.r))]
  
  # Caso directo
  m[, case := fifelse(set_usr == set_ref, "match", NA_character_)]
  
  # Complemento
  m[is.na(case), `:=`(
    A1c = .comp_map[A1.s],
    A2c = .comp_map[A2.s],
    set_usr_c = paste(pmin(.comp_map[A1.s], .comp_map[A2.s]),
                      pmax(.comp_map[A1.s], .comp_map[A2.s]))
  )]
  m[is.na(case) & set_usr_c == set_ref, case := "comp"]
  
  # Palindròmics
  m[, pal := .is_palindromic(A1.s, A2.s)]
  
  # Mantener (match o comp) y palindrómicos según opción
  keep <- m[case %chin% c("match","comp")]
  if (drop_palindromic) keep <- keep[pal == FALSE]
  
  # Aplicar complemento cuando haga falta
  keep[, was_comp := FALSE]
  keep[case == "comp", `:=`(
    A1.s = .comp_map[A1.s],
    A2.s = .comp_map[A2.s],
    was_comp = TRUE
  )]
  
  # Swap si orientaciones invertidas respecto a ref
  need_swap <- (keep$A1.s == keep$A2.r & keep$A2.s == keep$A1.r)
  keep[, was_swap := FALSE]
  if (any(need_swap)) {
    tmp <- keep$A1.s[need_swap]
    keep$A1.s[need_swap] <- keep$A2.s[need_swap]
    keep$A2.s[need_swap] <- tmp
    keep$was_swap[need_swap] <- TRUE
  }
  
  # ----- ORDEN Y SALIDA SEGURA -----
  
  # Orden de salida = orden del usuario (idx)
  setorder(keep, idx)
  
  # Construye .bim armonizado
  out <- keep[, .(CHR, SNP, CM, BP, A1 = A1.s, A2 = A2.s)]
  if (!nrow(out)) stop("❌ Harmonization fails: cap SNP vàlid")
  
  # Determina ruta final SIN duplicar extensión y escribe
  .out_file <- if (grepl("\\.bim$", out_bim_path, ignore.case = TRUE)) {
    out_bim_path
  } else {
    paste0(out_bim_path, ".bim")
  }
  fwrite(out, .out_file, sep = "\t", col.names = FALSE)
  
  # Comprueba existencia
  if (!file.exists(.out_file)) {
    stop("❌ No se generó el fichero: ", normalizePath(.out_file, mustWork = FALSE))
  }
  

  # Motivos de descarte (para informe)
  m[, remove_reason := NA_character_]
  m[!(case %chin% c("match","comp")), remove_reason := "no_match"]
  if (drop_palindromic) m[pal == TRUE, remove_reason := "palindromic"]
  
  not_in_ref <- usr[!SNP %in% ref$SNP]
  not_in_ref[, `:=`(remove_reason = "not_in_reference",
                    A1.s = A1.s, A2.s = A2.s)]
  
  removed <- rbind(
    m[!SNP %in% keep$SNP, .(CHR, SNP, CM, BP, A1 = A1.s, A2 = A2.s, remove_reason)],
    not_in_ref[, .(CHR, SNP, CM, BP, A1 = A1.s, A2 = A2.s, remove_reason)],
    fill = TRUE
  )
  
  # Resumen
  result <- list(
    n_input    = nrow(usr),
    n_kept     = nrow(out),
    n_removed  = nrow(removed),
    harmonized = out,
    removed    = removed,
    stats      = list(
      palindromic_dropped = sum(m$pal, na.rm = TRUE) * as.integer(drop_palindromic),
      matched             = sum(m$case == "match", na.rm = TRUE),
      complemented        = sum(m$case == "comp",  na.rm = TRUE)
    )
  )
  
  if (interactive()) {
    message(sprintf("✅ Harmonization completed: %d/%d SNPs retained",
                    result$n_kept, result$n_input))
  }
  
  invisible(result)
}
