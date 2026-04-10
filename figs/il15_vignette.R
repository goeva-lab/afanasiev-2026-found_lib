suppressPackageStartupMessages({
  library(glmGamPoi)
  library(data.table)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
})

cache <- function(x, path) {
  if (file.exists(path)) {
    readRDS(path)
  } else {
    saveRDS(x, file = path)
    x
  }
}

sce <- local({
  sce <- cbind(
    ... = purrr::map(c("pbs", "il15"), function(e) {
      anndataR::read_h5ad(sprintf("./.cache/vignette/%s.h5ad"), as = "SingleCellExperiment")
    }),
  )
  assayNames(sce) <- "counts"
  sce$cytokine <- relevel(sce$cytokine, ref = "PBS")
  outs <- found::HiDDENg(
    sce,
    k = 50L,
    group_by = "cell_type",
    cond_col = "cytokine",
    control_val = "PBS",
    X = t(assay(sce))
  )
  colData(sce)[["HiDDEN_p_hat"]] <- outs[["p_hat"]]
  colData(sce)[["HiDDEN_labs"]] <- outs[["labs"]]
  sce
}) |>
  cache(".cache/vignette/il15.SingleCellExperiment.rds")


de.res.bin <- list(
  "w/o HiDDEN" = sce,
  "after HiDDEN" = sce[, sce$HiDDEN_labs == sce$cytokine]
) |>
  purrr::map(function(sce) {
    sce$cell_type |>
      levels() |>
      setNames(nm = _) |>
      purrr::map(
        function(ct, fit.full) {
          ct_mask <- fit.full[["data"]]$cell_type == ct

          tryCatch(
            test_de(
              glm_gp(
                fit.full[["data"]][, ct_mask],
                design = ~cytokine,
                size_factors = fit.full[["size_factors"]][ct_mask],
                overdispersion = fit.full[["overdispersions"]]
              ),
              cond(cytokine = "IL-15") - cond(cytokine = "PBS"),
              compute_lfc_se = TRUE,
              max_lfc = Inf
            ),
            error = function(e) data.table()
          )
        },
        fit.full = glm_gp(
          pseudobulk(sce, vars(donor, cytokine, cell_type)),
          ~ cytokine + cell_type,
          size_factors = "normed_sum",
          verbose = TRUE
        )
      ) |>
      rbindlist(idcol = "cell_type")
  }) |>
  rbindlist(idcol = "method") |>
  cache(".cache/vignette/de_bin.data.table.rds")

de.res.cnt <- sce$cell_type |>
  levels() |>
  setNames(nm = _) |>
  purrr::map(
    function(ct, fit.full) {
      ct_mask <- sce$cell_type == ct
      fit.ct <- glm_gp(
        fit.full[["data"]][, ct_mask],
        design = fit.full[["model_matrix"]][ct_mask, ],
        size_factors = fit.full[["size_factors"]][ct_mask],
        overdispersion = fit.full[["overdispersions"]]
      )

      tryCatch(
        test_de(fit.ct, HiDDEN_p_hat, compute_lfc_se = TRUE, max_lfc = Inf),
        error = function(e) data.table()
      )
    },
    fit.full = glm_gp(sce, ~HiDDEN_p_hat, verbose = TRUE)
  ) |>
  rbindlist(idcol = "cell_type") |>
  cache(".cache/vignette/de_reg.data.table.rds")


de.res.cnt[, let(
  cell_type = as.factor(cell_type)
)]
de.res.bin[, let(
  cell_type = as.factor(cell_type),
  method = as.factor(method)
)]


ct.oi.cnt <- c("NK", "CD8 Memory", "CD14 Mono")
ct.oi.bin <- c("B Naive", "CD4 Memory", "MAIT")

vign.pl <- wrap_plots(
  as.data.table(colData(sce))[cell_type %in% ct.oi.cnt] |>
    ggplot() +
    aes(
      x = cytokine,
      y = HiDDEN_p_hat,
    ) +
    geom_violin(
      position = position_identity(),
      fill = "grey",
      scale = "count",
    ) +
    ggh4x::facet_nested(
      cols = vars(forcats::fct_relevel(cell_type, ct.oi.cnt), cytokine),
      scales = "free"
    ) +
    labs(fill = NULL, y = "HiDDEN scores", x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    guides(fill = guide_legend(position = "bottom")),
  de.res.cnt[cell_type %in% ct.oi.cnt] |>
    ggplot() +
    aes(
      x = pmax(pmin(lfc, 10), -10),
      y = -log10(adj_pval),
      colour = fcase(
        adj_pval >= 0.05 , "not significant" , # nolint: commas_linter.
        lfc > 0          , "up"              , # nolint: commas_linter.
        default = "down"
      ) |>
        forcats::fct_relevel("not significant")
    ) +
    geom_point(shape = ".") +
    geom_label(
      aes(
        label = n,
        x = fifelse(lfc, 7.5, -7.5),
        colour = fifelse(lfc, "up", "down")
      ),
      de.res.cnt[
        (cell_type %in% ct.oi.cnt) & (adj_pval < 0.05),
        .(n = .N),
        by = .(lfc > 0, cell_type)
      ],
      y = 150,
      show.legend = FALSE
    ) +
    facet_grid(cols = vars(forcats::fct_relevel(cell_type, ct.oi.cnt))) +
    labs(
      colour = NULL,
      x = "change coefficient w/ HiDDEN score",
      y = expression(-log[10](pval))
    ) +
    guides(colour = guide_legend(override.aes = list(shape = 16))) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank()
    ),
  dcast(
    rbind(
      "continuous" = de.res.cnt[
        (adj_pval < 0.05) & (cell_type %in% ct.oi.cnt)
      ],
      "classic" = de.res.bin[
        (adj_pval < 0.05) &
          (cell_type %in% ct.oi.cnt) &
          (method == "w/o HiDDEN")
      ][, -c("method")],
      idcol = "method"
    ),
    name + cell_type ~ method,
    value.var = c("adj_pval", "lfc")
  )[,
    .(n = .N),
    by = .(
      cell_type,
      annot_pval = forcats::fct_relevel(
        fcase(
          `adj_pval_continuous` < 0.05 & `adj_pval_classic` < 0.05 , # nolint: commas_linter.
          "significant in\nboth analyses"                          , # nolint: commas_linter.
          `adj_pval_classic` < 0.05                                , # nolint: commas_linter.
          "not significant\nw/ regression"                         , # nolint: commas_linter.
          `adj_pval_continuous` < 0.05                             , # nolint: commas_linter.
          "only significant\nw/ regression"                        , # nolint: commas_linter.
          default = "not significant"
        ),
        "not significant\nw/ regression",
        after = Inf
      )
    )
  ] |>
    dcast(cell_type ~ annot_pval, value.var = "n", fill = 0) |>
    melt(
      id.vars = "cell_type",
      value.name = "n",
      variable.name = "annot_pval"
    ) |>
    ggplot() +
    aes(
      x = cell_type,
      fill = annot_pval,
      y = n
    ) +
    geom_col(position = position_dodge(), linewidth = 0) +
    scale_y_reverse() +
    facet_wrap(
      vars(forcats::fct_relevel(cell_type, ct.oi.bin)),
      nrow = 1,
      scales = "free_x"
    ) +
    labs(x = NULL, y = "DEG count", fill = NULL) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_fill_brewer(palette = "Dark2"),
  as.data.table(colData(sce))[cell_type %in% ct.oi.bin] |>
    ggplot() +
    aes(
      x = cytokine,
      y = HiDDEN_p_hat,
      fill = fifelse(
        (cytokine == "PBS") | (HiDDEN_labs != cytokine),
        "HiDDEN:\nunaffected",
        "HiDDEN:\naffected"
      )
    ) +
    geom_violin(
      position = position_identity(),
      alpha = 0.8,
      scale = "count",
    ) +
    ggh4x::facet_nested(
      cols = vars(forcats::fct_relevel(cell_type, ct.oi.bin), cytokine),
      scales = "free"
    ) +
    labs(fill = NULL, y = "HiDDEN scores", x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  de.res.bin[cell_type %in% ct.oi.bin] |>
    ggplot() +
    aes(
      x = pmax(pmin(lfc, 10), -10),
      y = -log10(adj_pval),
      colour = fcase(
        adj_pval >= 0.05 , "not significant" , # nolint: commas_linter.
        lfc > 0          , "up"              , # nolint: commas_linter.
        default = "down"
      ) |>
        forcats::fct_reorder(adj_pval, .desc = TRUE)
    ) +
    geom_point(shape = ".") +
    geom_label(
      aes(
        label = n,
        x = fifelse(lfc, 7.5, -7.5),
        colour = fifelse(lfc, "up", "down")
      ),
      de.res.bin[
        (cell_type %in% ct.oi.bin) & adj_pval < 0.05,
        .(n = .N),
        by = .(lfc > 0, cell_type, method)
      ],
      y = 5,
      show.legend = FALSE
    ) +
    facet_grid(
      cols = vars(forcats::fct_relevel(cell_type, ct.oi.bin)),
      rows = vars(forcats::fct_relevel(method, "w/o HiDDEN")),
      scales = "free_y"
    ) +
    labs(
      colour = NULL,
      x = expression(log[2] ~ fold ~ change, ""),
      y = expression(-log[10](pval))
    ) +
    lims(x = c(-10, 10)) +
    guides(colour = guide_legend(override.aes = list(shape = 16))) +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = element_blank()
    ),
  dcast(
    de.res.bin[adj_pval < 0.05 & cell_type %in% ct.oi.bin],
    name + cell_type ~ method,
    value.var = c("adj_pval", "lfc")
  )[,
    .(n = .N),
    by = .(
      cell_type,
      annot_pval = forcats::fct_relevel(
        fcase(
          `adj_pval_after HiDDEN` < 0.05 & `adj_pval_w/o HiDDEN` < 0.05 , # nolint: commas_linter.
          "significant in\nboth analyses"                               , # nolint: commas_linter.
          `adj_pval_w/o HiDDEN` < 0.05                                  , # nolint: commas_linter.
          "not significant\nafter HiDDEN"                               , # nolint: commas_linter.
          `adj_pval_after HiDDEN` < 0.05                                , # nolint: commas_linter.
          "only significant\nafter HiDDEN"                              , # nolint: commas_linter.
          default = "not significant"
        ),
        "not significant\nafter HiDDEN",
        after = Inf
      )
    )
  ] |>
    dcast(cell_type ~ annot_pval, value.var = "n", fill = 0) |>
    melt(
      id.vars = "cell_type",
      value.name = "n",
      variable.name = "annot_pval"
    ) |>
    ggplot() +
    aes(
      x = cell_type,
      fill = annot_pval,
      y = n
    ) +
    geom_col(position = position_dodge(), linewidth = 0) +
    scale_y_reverse() +
    facet_wrap(
      vars(forcats::fct_relevel(cell_type, ct.oi.bin)),
      nrow = 1,
      scales = "free_x"
    ) +
    labs(x = NULL, y = "DEG count", fill = NULL) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_fill_brewer(palette = "Dark2"),
  heights = c(1, 1, 0.6, 1, 2, 0.6),
  design = "A\nB\nC\nD\nE\nF"
) &
  theme(legend.key.spacing.y = unit(0.1, "lines"))


if (!interactive()) {
  ggsave("fig_1d.svg", vign.pl, h = 3465, w = 2325, units = "px")
}
