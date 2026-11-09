# Smoothers in the Field: GAMs vs ASReml for spatial analysis of field trials
# Companion script for index.Rmd (single environment, LOC1).

library(mgcv)
library(asreml)      # ASReml-R (requires a license)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(viridis)


# Load and prepare data ---------------------------------------------------

yield_data <- read.csv("FARGO_pREP_2025-04-18.csv")

# Sorted by ROW then COLUMN; AR_start marks the first plot of each row so
# the AR(1) residual process restarts at every row in bam().
yield_data_single_gam <- yield_data |>
	filter(LOCATION == "LOC1") |>
	mutate(
		TREATMENT  = as.factor(TREATMENT),
		ROW_fc     = as.factor(ROW),
		COLUMN_fc  = as.factor(COLUMN)
	) |>
	arrange(ROW, COLUMN) |>
	group_by(ROW) |>
	mutate(AR_start = row_number() == 1) |>
	ungroup()

# ASReml needs the data sorted to match ar1(ROW_fc):ar1(COLUMN_fc).
# This is the same plot order as yield_data_single_gam.
yield_data_single_asreml <- yield_data_single_gam |>
	arrange(ROW_fc, COLUMN_fc)


# Helper functions --------------------------------------------------------

# Add confidence intervals
add_confidence_intervals <- function(tbl, pred_col, se_col, level = 0.95) {
	z <- qnorm(1 - (1 - level) / 2)
	tbl |>
		mutate(
			lower = {{ pred_col }} - z * {{ se_col }},
			upper = {{ pred_col }} + z * {{ se_col }}
		)
}

# Caterpillar plot of predicted treatment values with CIs
plot_treatment_caterpillar <- function(tbl, reorder = TRUE) {
	tbl <- if (reorder) {
		tbl |> mutate(TREATMENT = reorder(TREATMENT, -predicted_value))
	} else tbl

	ggplot(tbl, aes(x = TREATMENT, y = predicted_value)) +
		geom_point(size = 2.5) +
		geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
		labs(
			title = "Treatment Effects (GAM, Adjusted for Spatial Trends)",
			x = "Treatment",
			y = "Predicted Yield ± CI"
		) +
		theme_minimal(base_size = 10) +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
			plot.title = element_text(face = "bold"),
			panel.grid.major.x = element_blank()
		)
}

# Histogram, residuals vs fitted, and QQ plot side by side.
# res_tbl must contain the columns `fitted` and `residual`.
plot_residual_diagnostics <- function(res_tbl, model_label) {
	hist_plot <- ggplot(res_tbl, aes(x = residual)) +
		geom_histogram(binwidth = 0.3, fill = "steelblue", color = "black") +
		labs(
			title = paste0("Histogram of Residuals (", model_label, ")"),
			x = "Residual",
			y = "Frequency"
		) +
		theme_minimal(base_size = 10)

	resid_fit_plot <- ggplot(res_tbl, aes(x = fitted, y = residual)) +
		geom_point(alpha = 0.8, color = "black") +
		geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
		labs(
			x = "Fitted Values",
			y = "Residuals",
			title = paste0("Residuals vs Fitted Values (", model_label, ")")
		) +
		theme_minimal(base_size = 10)

	qq_plot <- ggplot(res_tbl, aes(sample = residual)) +
		stat_qq(color = "darkgreen") +
		stat_qq_line(color = "black", linewidth = 0.8) +
		labs(
			title = paste0("QQ Plot of Residuals (", model_label, ")"),
			x = "Theoretical Quantiles",
			y = "Sample Quantiles"
		) +
		theme_minimal(base_size = 10)

	hist_plot | resid_fit_plot | qq_plot
}

# Residuals on the field layout (purple = negative, yellow = positive)
plot_residual_heatmap <- function(res_tbl, model_label) {
	ggplot(res_tbl, aes(x = COLUMN, y = ROW, fill = residual)) +
		geom_tile(color = "white") +
		scale_fill_viridis(
			option = "viridis",  # purple → green → yellow
			name = "Residual",
			direction = 1
		) +
		scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
		scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
		labs(
			title = paste0("Spatial Heatmap of Residuals (", model_label, ")"),
			x = "Column",
			y = "Row"
		) +
		theme_minimal(base_size = 13) +
		theme(
			plot.title = element_text(face = "bold", hjust = 0.5),
			legend.position = "right"
		)
}

# Top-N agreement between two rankings: % of genotypes selected by both
# models and Kendall's tau among those shared genotypes.
compare_top_n <- function(rank_tbl, top_n) {
	top_gam    <- rank_tbl$TREATMENT[rank_tbl$rank_gam_total <= top_n]
	top_asreml <- rank_tbl$TREATMENT[rank_tbl$rank_asreml <= top_n]
	shared     <- intersect(top_gam, top_asreml)
	shared_tbl <- dplyr::filter(rank_tbl, TREATMENT %in% shared)

	data.frame(
		Top = top_n,
		kendall_tau = round(cor(shared_tbl$rank_gam_total, shared_tbl$rank_asreml, method = "kendall"), 3),
		shared_pct = round(100 * length(shared) / top_n, 1)
	)
}


# Estimate rho ------------------------------------------------------------

# bam() does not estimate the AR(1) parameter, so we pick the value with the
# lowest REML score. This fits 18 models and takes a minute or two.
test_rho <- function(rho_val) {
	bam(
		YIELD ~
			s(ROW, bs = "ps") + s(COLUMN, bs = "ps") +
			s(ROW_fc, bs = "re") + s(COLUMN_fc, bs = "re") +
			s(TREATMENT, bs = "re") +
			te(ROW_fc, COLUMN_fc, bs = c("re", "re")),
		rho = rho_val,
		AR.start = yield_data_single_gam$AR_start,
		data = yield_data_single_gam,
		method = "fREML"
	)$gcv.ubre
}

rhos <- seq(0.1, 0.95, by = 0.05)
scores <- sapply(rhos, test_rho)

rho_search <- data.frame(rho = rhos, reml_score = scores) |>
	arrange(reml_score)

head(rho_search, 5)
best_rho <- rho_search$rho[1]  # 0.55 for LOC1


# Fit GAM model -----------------------------------------------------------

gam_m2 <- bam(
	YIELD ~
		s(ROW, bs = "ps") + s(COLUMN, bs = "ps") +
		s(ROW_fc, bs = "re") + s(COLUMN_fc, bs = "re") +
		s(TREATMENT, bs = "re") +
		te(ROW_fc, COLUMN_fc, bs = c("re", "re")),
	rho = best_rho,
	AR.start = yield_data_single_gam$AR_start,
	data = yield_data_single_gam,
	method = "fREML"
)


# GAM genotype predictions ------------------------------------------------

# Smooth trends are evaluated at the field center, and the random row,
# column, and row x column effects are excluded (set to zero). ROW_fc and
# COLUMN_fc must be present in newdata, so any valid level is used.
newdata_gam <- tibble::tibble(
	TREATMENT = levels(yield_data_single_gam$TREATMENT),
	ROW = mean(yield_data_single_gam$ROW),
	COLUMN = mean(yield_data_single_gam$COLUMN),
	ROW_fc = yield_data_single_gam$ROW_fc[1],
	COLUMN_fc = yield_data_single_gam$COLUMN_fc[1]
)

preds_treatment <- predict(
	gam_m2,
	newdata = newdata_gam,
	se.fit = TRUE,
	type = "response",
	exclude = c("s(ROW_fc)", "s(COLUMN_fc)", "te(ROW_fc,COLUMN_fc)")
)

BLUPs_gam_single <- newdata_gam |>
	mutate(
		predicted_value = as.numeric(preds_treatment$fit),
		se = as.numeric(preds_treatment$se.fit)
	) |>
	add_confidence_intervals(pred_col = predicted_value, se_col = se) |>
	select(TREATMENT, predicted_value, se, lower, upper)

BLUPs_gam_single |> head(10)

plot_treatment_caterpillar(BLUPs_gam_single)


# GAM residual diagnostics ------------------------------------------------

yield_data_gam_res <- yield_data_single_gam |>
	mutate(
		fitted   = as.numeric(predict(gam_m2, type = "response")),
		residual = as.numeric(residuals(gam_m2, type = "response"))
	)

plot_residual_diagnostics(yield_data_gam_res, "GAM")


# Fit ASReml model --------------------------------------------------------

asreml_m2 <- asreml(
	fixed = YIELD ~ ROW + COLUMN,
	random = ~ spl(ROW) + spl(COLUMN) + TREATMENT,
	residual = ~ ar1(ROW_fc):ar1(COLUMN_fc),
	data = yield_data_single_asreml,
	trace = FALSE
)

# Refit from the current estimates until the model converges
for (i in 1:5) {
	if (isTRUE(asreml_m2$converge)) break
	asreml_m2 <- update(asreml_m2)
}
if (!isTRUE(asreml_m2$converge)) warning("ASReml model did not converge")


# ASReml genotype predictions ---------------------------------------------

pred_asreml <- predict(asreml_m2, classify = "TREATMENT", sed = TRUE)$pvals

pred_asreml |>
	add_confidence_intervals(pred_col = predicted.value, se_col = std.error) |>
	head(10)


# ASReml residual diagnostics ---------------------------------------------

yield_data_single_asreml_res <- yield_data_single_asreml |>
	mutate(
		fitted = fitted(asreml_m2),
		residual = residuals(asreml_m2, type = "response")
	)

plot_residual_diagnostics(yield_data_single_asreml_res, "ASReml")


# Residual heatmaps -------------------------------------------------------

(plot_residual_heatmap(yield_data_gam_res, "GAM") |
 	plot_residual_heatmap(yield_data_single_asreml_res, "ASReml")) +
	plot_layout(guides = "collect") &
	theme(legend.position = "right")


# Correlation of predicted values -----------------------------------------

# Match predictions by genotype instead of relying on row order
pred_comparison <- BLUPs_gam_single |>
	select(TREATMENT, pred_gam = predicted_value, se_gam = se) |>
	inner_join(
		pred_asreml |>
			mutate(TREATMENT = as.character(TREATMENT)) |>
			select(TREATMENT, pred_asreml = predicted.value, se_asreml = std.error),
		by = "TREATMENT"
	)

stopifnot(nrow(pred_comparison) == nlevels(yield_data_single_gam$TREATMENT))

cor(pred_comparison$pred_gam, pred_comparison$pred_asreml)


# Agreement in genotype rankings ------------------------------------------

# rank 1 = highest predicted yield
rank_comparison <- pred_comparison |>
	mutate(
		rank_gam_total = rank(-pred_gam),
		rank_asreml = rank(-pred_asreml),
		diff = rank_gam_total - rank_asreml
	)

# Largest rank differences
rank_comparison |> arrange(desc(abs(diff))) |> head(10)

# Kendall's tau over all genotypes
cor(rank_comparison$rank_gam_total, rank_comparison$rank_asreml, method = "kendall")

# Top-N selection agreement
rankings <- do.call(
	rbind,
	lapply(c(10, 20, 30, 40, 50), compare_top_n, rank_tbl = rank_comparison)
)

rankings |> kable()

ggplot(rank_comparison, aes(x = rank_asreml, y = rank_gam_total)) +
	geom_point(size = 2.5, color = "darkgreen") +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
	labs(
		title = "Rank Comparison: ASReml vs GAM (Predicted Genotype Values)",
		x = "ASReml Rank",
		y = "GAM Rank"
	) +
	theme_minimal(base_size = 13)


# Standard errors of predictions ------------------------------------------

se_long <- pred_comparison |>
	select(TREATMENT, se_gam, se_asreml) |>
	pivot_longer(
		cols = c(se_gam, se_asreml),
		names_to = "model",
		values_to = "se"
	) |>
	mutate(
		model = recode(model, se_gam = "GAM", se_asreml = "ASReml"),
		TREATMENT = factor(TREATMENT, levels = unique(TREATMENT))  # preserve order
	)

ggplot(se_long, aes(x = TREATMENT, y = se, group = model, color = model)) +
	geom_line(linewidth = 1) +
	geom_point(size = 1.5) +
	scale_color_manual(values = c("GAM" = "#1b9e77", "ASReml" = "#d95f02")) +
	labs(
		title = "Standard Errors of Treatment Predictions by Model",
		x = "Genotype (Treatment)",
		y = "Standard Error",
		color = "Model"
	) +
	theme_minimal(base_size = 13) +
	theme(
		axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
		plot.title = element_text(face = "bold"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)
