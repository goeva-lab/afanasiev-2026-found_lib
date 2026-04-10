suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

gst <- function(dir, filter_regex = ".*") {
  out.l <- nc::capture_first_vec(
    stringr::str_subset(
      list.files(
        file.path("cache", dir),
        pattern = "out.tsv",
        recursive = TRUE
      ),
      stringr::regex(filter_regex)
    ),
    emb = "[^/]+",
    "/",
    k = "k<?=[0-9]+",
    "/",
    reg = "(?|(?|logit|svm|rf)_[^/]+)?",
    "/?",
    bin = "(?|gmm|kmeans)?",
    "/?",
    grp = "(?|[^/]+)?",
    "/?",
    "out.tsv"
  ) |>
    split(by = c("emb", "k", "reg", "bin"), keep.by = FALSE, flatten = FALSE) |>
    purrr::imap(function(l, emb) {
      purrr::imap(l, function(l, k) {
        purrr::imap(l, function(l, reg) {
          bin.l <- purrr::imap(l, function(grps, bin) {
            idv <- c(emb, k, reg, bin)
            idv <- idv[nzchar(idv)]
            grp.dt <- rbindlist(
              purrr::map(
                grps[, grp[nzchar(grp)]],
                function(e) {
                  fread(stringr::str_c(
                    c("cache", dir, idv, e, "out.tsv"),
                    collapse = .Platform[["file.sep"]]
                  ))
                }
              )
            )
            list(
              "all" = setnames(
                fread(stringr::str_c(
                  c("cache", dir, idv, "out.tsv"),
                  collapse = .Platform[["file.sep"]]
                ))[grp.dt[, cell], on = .(cell)],
                colnames(grp.dt[, -c("cell")]),
                stringr::str_c(c("all", idv), collapse = "/")
              ),
              "grp" = setnames(
                grp.dt,
                colnames(grp.dt[, -c("cell")]),
                stringr::str_c(c("grp", idv), collapse = "/")
              )
            )
          })
          if ((length(bin.l) == 1) && !nzchar(names(bin.l)[[1]])) {
            bin.l[[1]]
          } else {
            bin.l
          }
        })
      })
    })

  while (class(out.l[[1]])[[1]] == "list") {
    out.l <- purrr::list_flatten(out.l)
  }

  purrr::reduce(out.l, function(lhs, rhs) lhs[rhs, on = .(cell)])
}

sdt <- function(dists.dt, col) {
  DT <- dists.dt[,
    transpose(as.data.table(
      stringr::str_split(get(col), stringr::fixed("/"))
    ))
  ]
  cbind(
    dists.dt,
    setnames(
      DT,
      sprintf(c("%s.grp", "%s.emb", "%s.k", "%s.reg", "%s.bin"), col)[
        seq_len(ncol(DT))
      ]
    )[,
      sprintf("%s.k", col) := as.integer(
        stringr::str_split_i(get(sprintf("%s.k", col)), "=", 2)
      )
    ][]
  )
}

emb.ord <- local({
  emb.ord.naive <- c(
    "log_pca",
    "sd_log_pca",
    "rs_nmf",
    "log_nmf",
    "sd_nmf",
    "sd_rs_nmf",
    "sd_log_nmf",
    "res_pca",
    "sd_res_pca"
  )
  c(emb.ord.naive, sprintf("%s_harmony", emb.ord.naive), "scVI")
})

cache <- function(x, path) {
  if (file.exists(path)) {
    readRDS(path)
  } else {
    saveRDS(x, file = path)
    x
  }
}


META.LIST <- local({
  ad <- reticulate::import("anndata", convert = FALSE)
  data_conf <- reticulate::import_from_path(
    "bench",
    convert = FALSE
  )$DATA_CONF
  purrr::map(
    setNames(
      nm = reticulate::py_to_r(
        reticulate::import_builtins()$list(data_conf$keys())
      )
    ),
    function(k) {
      entry <- data_conf[k]
      adata <- ad$read_h5ad(
        sprintf("cache/data/%s.h5ad", reticulate::py_to_r(entry[0])),
        backed = "r"
      )
      meta <- adata$obs$iloc[entry[1](adata$obs)]
      m.dt <- data.table(
        cell = reticulate::py_to_r(meta$index$to_numpy()),
        cond = fifelse(
          reticulate::py_to_r(entry[2](meta)),
          "case",
          "control"
        ),
        type = reticulate::py_to_r(entry[3](meta)$to_numpy()),
        batch = reticulate::py_to_r(entry[4](meta)$to_numpy())
      )
      attr(m.dt, "n.var") <- reticulate::py_to_r(adata$n_vars)
      m.dt
    }
  )
}) |>
  cache(".cache/meta.list.rds")

reg.dt.l <- parallel::mclapply(
  setNames(nm = names(META.LIST)),
  function(ds) gst(file.path("reg", ds)),
  mc.cores = parallel::detectCores()
) |>
  cache(".cache/bench/reg.list.rds")

bin.dt.l <- parallel::mclapply(
  setNames(nm = names(META.LIST)),
  function(ds) gst(file.path("bin", ds)),
  mc.cores = parallel::detectCores()
) |>
  cache(".cache/bench/bin.list.rds")

score.null.dt <- rbindlist(
  purrr::pmap(
    list(reg.dt.l, META.LIST),
    function(reg.dt, meta.dt) {
      scores <- stringr::str_subset(
        colnames(reg.dt)[-1],
        stringr::fixed("_shuf"),
        negate = TRUE
      ) |>
        setNames(nm = _) |>
        parallel::mclapply(
          function(e, reg.dt) {
            approxOT::wasserstein(
              reg.dt[[e]],
              reg.dt[[sprintf("%s_shuf", e)]],
              method = "univariate",
              p = 1
            )
          },
          mc.cores = parallel::detectCores(),
          reg.dt = reg.dt
        ) |>
        c(recursive = TRUE)

      sdt(data.table(meth = names(scores), score = scores), "meth")
    }
  ),
  idcol = "dataset"
) |>
  cache(".cache/score_null.data.table.rds")

score.dist.dt <- rbindlist(
  purrr::pmap(
    list(reg.dt.l, bin.dt.l, META.LIST),
    function(reg.dt, bin.dt, meta.dt) {
      scores <- colnames(bin.dt)[-1] |>
        setNames(nm = _) |>
        parallel::mclapply(
          function(e, reg.dt, bin.dt, cond_mask) {
            regs <- split(
              reg.dt[[stringr::str_extract(e, "(.+)/[^/]+$", 1)]],
              cond_mask
            )
            bins <- bin.dt[[e]][cond_mask]
            regs.case <- split(regs[["TRUE"]], bins)
            if (length(regs.case) == 1) {
              NA # early exit if no or total relabeling
            } else {
              approxOT::wasserstein(
                regs.case[["TRUE"]],
                regs.case[["FALSE"]],
                p = 1,
                method = "univariate"
              ) -
                approxOT::wasserstein(
                  regs[["TRUE"]],
                  regs.case[["FALSE"]],
                  p = 1,
                  method = "univariate"
                )
            }
          },
          mc.cores = parallel::detectCores(),
          reg.dt = reg.dt,
          # ensure conserved order across data.tables
          # should be unnecessary? but good sanity check
          bin.dt = bin.dt[reg.dt[, cell], on = .(cell)],
          cond_mask = meta.dt[
            reg.dt[, cell],
            cond != "control",
            on = .(cell)
          ]
        ) |>
        c(recursive = TRUE)

      sdt(data.table(meth = names(scores), score = scores), "meth")
    }
  ),
  idcol = "dataset"
) |>
  cache(".cache/score_dist.data.table.rds")


emb.diff.pl <- score.null.dt[meth.reg == "logit_lbfgs_nol1"][
  order(score),
  last(.SD),
  by = .(
    dataset,
    meth.emb = stringr::str_remove(meth.emb, stringr::fixed("_harmony")),
    meth.grp
  )
] |>
  ggplot() +
  aes(
    x = score,
    y = forcats::fct(
      dataset,
      names(reg.dt.l)[order(
        purrr::map_vec(reg.dt.l, nrow),
        decreasing = TRUE
      )]
    ),
    shape = meth.grp
  ) +
  geom_violin(
    aes(fill = meth.grp, colour = after_scale(fill)),
    position = position_identity(),
    linewidth = 0.5,
    alpha = 0.7,
    scale = "width"
  ) +
  geom_point(aes(colour = forcats::fct(meth.emb, emb.ord))) +
  lims(x = c(0, NA)) +
  labs(y = NULL, colour = NULL, shape = NULL, fill = NULL)
k.diff.pl <- score.dist.dt[meth.reg == "logit_lbfgs_nol1" & !is.na(score)][
  order(meth.k),
  .(norm.score = score / score[[1]], meth.k),
  by = .(dataset, meth.emb, meth.grp, meth.bin)
] |>
  ggplot() +
  aes(
    x = meth.k,
    y = norm.score,
    colour = forcats::fct(
      dataset,
      names(reg.dt.l)[order(purrr::map_vec(reg.dt.l, nrow))]
    ),
    fill = after_scale(colour),
  ) +
  geom_smooth(formula = y ~ log10(x), method = "loess") +
  facet_grid(cols = vars(meth.grp), scales = "free_y") +
  geom_hline(
    yintercept = 1,
    alpha = 0.5,
    colour = "red",
    linetype = "dashed"
  ) +
  scale_y_log10() +
  scale_x_log10(breaks = c(2:5, 10, 30, 50)) +
  scale_colour_viridis_d(option = "turbo") +
  labs(
    linetype = NULL,
    colour = NULL,
    x = "k (log-scale)",
    y = "score fold change (log-scale)"
  ) +
  theme(legend.key.spacing.y = unit(0.1, "lines"))


reg.cmp.pl <- purrr::pmap(
  list(
    split(
      score.dist.dt[
        meth.bin == "kmeans" &
          stringr::str_detect(meth.emb, stringr::fixed("sd_log_pca"))
      ][
        order(score, na.last = FALSE),
        last(.SD),
        by = .(dataset, meth.reg)
      ],
      by = "dataset"
    )[names(META.LIST)],
    reg.dt.l,
    META.LIST,
    bin.dt.l
  ),
  function(col.dt, reg.dt, meta.dt, bin.dt) {
    bin.cols <- col.dt[, meth]
    reg.cols <- stringr::str_sub(bin.cols, end = -8)
    melt(
      meta.dt[
        reg.dt[, mget(c("cell", reg.cols))],
        on = .(cell)
      ][
        bin.dt[, mget(c("cell", bin.cols))],
        on = .(cell)
      ],
      id.vars = colnames(meta.dt),
      measure.vars = measure(
        meth,
        value.name = function(e) fifelse(nzchar(e), "kmeans.bin", "reg.out"),
        pattern = "([^/]+/[^/]+/[^/]+/[^/]+)(/kmeans)?"
      ),
    )
  }
) |>
  rbindlist(idcol = "dataset") |>
  ggplot() +
  aes(
    y = forcats::fct_relevel(cond, "control"),
    x = reg.out,
    fill = fifelse(
      kmeans.bin,
      "HiDDEN - affected",
      "HiDDEN - unaffected",
    ),
    colour = after_scale(fill)
  ) +
  geom_violin(
    scale = "count",
    linewidth = 0.2,
    alpha = 0.8,
    position = position_identity()
  ) +
  ggh4x::facet_grid2(
    rows = vars(forcats::fct(
      dataset,
      names(reg.dt.l)[order(purrr::map_vec(reg.dt.l, nrow))]
    )),
    cols = vars(stringr::str_extract(meth, "svm|logit|rf")),
    scales = "free_x",
    independent = "x",
  ) +
  scale_y_discrete(limits = rev) +
  ggh4x::scale_x_facet(COL < 3, limits = c(0, 1)) +
  labs(y = NULL, x = "HiDDEN score", fill = NULL) +
  guides(fill = guide_legend(position = "bottom")) +
  theme(legend.key.spacing.y = unit(0.1, "lines"))


bin.diff.pl <- dcast(
  score.dist.dt[meth.reg == "logit_lbfgs_nol1"],
  value.var = "score",
  dataset + meth.emb + meth.grp + meth.k ~ meth.bin
) |>
  na.omit() |>
  ggplot() +
  aes(
    x = kmeans - gmm,
    y = forcats::fct(
      stringr::str_remove(meth.emb, stringr::fixed("_harmony")),
      rev(emb.ord)
    ),
    colour = as.factor(meth.k),
    shape = meth.grp,
    alpha = stringr::str_ends(meth.emb, stringr::fixed("_harmony"))
  ) +
  geom_point() +
  scale_alpha_manual(
    labels = c("unadjusted", "harmonized"),
    values = c(0.8, 0.4)
  ) +
  facet_grid(
    rows = vars(forcats::fct(
      dataset,
      names(reg.dt.l)[order(purrr::map_vec(reg.dt.l, nrow))]
    ))
  ) +
  labs(
    colour = "k",
    shape = "grouping\nmode",
    alpha = "adjustment",
    y = NULL
  ) +
  scale_color_viridis_d(option = "turbo") +
  lims(x = c(-0.5, 0.5)) +
  theme(legend.key.spacing.y = unit(0.1, "lines"))


n.table <- rbind(
  rbindlist(
    META.LIST,
    idcol = "dataset"
  )[, .(n.obs = .N), by = .(grp = type, dataset)],
  data.table(
    dataset = names(META.LIST),
    n.obs = purrr::map_vec(META.LIST, nrow),
    grp = NA
  )
)[
  data.table(
    dataset = names(META.LIST),
    n.var = purrr::map_vec(META.LIST, attr, which = "n.var")
  ),
  on = .(dataset)
]
perf.dt.l <- purrr::map(
  setNames(nm = c("emb", "reg", "bin")),
  function(e) {
    n.table[
      fread(sprintf(".cache/perf/%s.tsv", e), na.strings = ""),
      on = .(dataset, grp)
    ][, let(
      total_time = lubridate::as.duration(total_time),
      peak_mem = as.numeric(peak_mem) / 1e6
    )][]
  }
)

perf.pl <- local({
  base.pl <- ggplot() +
    aes(
      fill = after_scale(colour),
      ymin = stage(after_stat = exp(as.numeric(ymin))),
      ymax = stage(after_stat = exp(as.numeric(ymax)))
    ) +
    geom_smooth(method = "gam", formula = log(y) ~ s(x), alpha = 0.2) +
    scale_x_log10() +
    labs(colour = NULL, fill = NULL, linetype = NULL) +
    theme(legend.key.spacing.y = unit(0.1, "lines"))
  emb.base <- add_gg(base.pl, perf.dt.l[["emb"]]) +
    aes(
      x = as.double(n.obs) * as.double(n.var),
      colour = forcats::fct(stringr::str_remove(emb, "_harmony"), emb.ord),
      linetype = forcats::fct_relevel(
        fifelse(endsWith(emb, "_harmony"), "harmonized", "unadjusted"),
        "unadjusted"
      )
    ) +
    labs(x = "size (obs * vars) (log-scale)") +
    scale_colour_brewer(palette = "Set3")
  reg.base <- add_gg(base.pl, perf.dt.l[["reg"]][!endsWith(reg, "_shuf")]) +
    aes(
      x = as.double(n.obs) * as.double(k),
      colour = stringr::str_extract(reg, stringr::regex("logit|rf|svm"))
    ) +
    labs(x = "size (obs * k) (log-scale)")
  bin.base <- add_gg(base.pl, perf.dt.l[["bin"]]) +
    aes(x = as.double(n.obs), colour = bin) +
    labs(x = "size (obs) (log-scale)")
  time.pl <- function(e) {
    e +
      aes(y = stage(total_time, after_stat = exp(as.numeric(y)))) +
      scale_y_time(limits = c(0, NA)) +
      labs(y = "time ([hh]:[mm]:[ss])") +
      guides(colour = "none", fill = "none", linetype = "none")
  }
  mem.pl <- function(e) e + aes(y = stage(peak_mem, after_stat = exp(y))) + labs(y = "peak memory (MBs)")

  wrap_plots(
    time.pl(emb.base),
    mem.pl(emb.base),
    time.pl(reg.base),
    mem.pl(reg.base),
    time.pl(bin.base),
    mem.pl(bin.base),
    design = "AB\nCD\nEF",
    heights = c(3, 2, 2)
  )
})


if (!interactive()) {
  ggsave(
    "fig_2.svg",
    wrap_plots(emb.diff.pl, k.diff.pl, design = "A\nB", heights = c(1, 1)),
    w = 3520,
    h = 2720,
    units = "px"
  )
  ggsave("fig_s1.svg", reg.cmp.pl, h = 5500, w = 4250, units = "px")
  ggsave("fig_s2.svg", bin.diff.pl, h = 3520, w = 2720, units = "px")
  ggsave("fig_s3.svg", perf.pl, h = 3520, w = 2720, units = "px")
}
