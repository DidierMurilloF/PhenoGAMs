library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)
library(kableExtra)

# Load and prepare data
yield_data <- read.csv("FARGO_pREP_2025-04-18.csv")

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


# Estimate Rho ------------------------------------------------------------

test_rho <- function(rho_val) {
	bam(
		YIELD ~ 
			s(ROW, bs = "ps") + s(COLUMN, bs = "ps") +
			s(ROW_fc, bs = "re") + s(COLUMN_fc, bs = "re") +
			s(TREATMENT, bs = "re"),
		data = yield_data_single_gam,
		method = "fREML",
		AR.start = yield_data_single_gam$AR_start,
		rho = rho_val
	)$gcv.ubre  
}

rhos <- seq(0.1, 0.95, by = 0.05)
scores <- sapply(rhos, test_rho)

data.frame(rho = rhos, score = scores) |> 
	arrange(score) |> 
	head(10)


# Fit GAM model
m2 <- bam(
	YIELD ~ 
		s(ROW, bs = "ps") + s(COLUMN, bs = "ps") +
		s(ROW_fc, bs = "re") + s(COLUMN_fc, bs = "re") +
		s(TREATMENT, bs = "re") +
		te(ROW_fc, COLUMN_fc, bs = c("re", "re")),
	rho = 0.55,
	AR.start = yield_data_single_gam$AR_start,
	data = yield_data_single_gam,
	method = "fREML"
)

# Center values for marginal prediction
mean_row <- mean(yield_data_single_gam$ROW)
mean_column <- mean(yield_data_single_gam$COLUMN)
row_fc_center <- as.factor(round(mean_row))
col_fc_center <- as.factor(round(mean_column))

# Generate prediction dataset
newdata_gam <- tibble::tibble(
	TREATMENT = levels(yield_data_single_gam$TREATMENT),
	ROW = mean_row,
	COLUMN = mean_column,
	ROW_fc = row_fc_center,
	COLUMN_fc = col_fc_center
)

# Predict treatment means at center of field
preds_treatment_centered <- predict(m2, newdata = newdata_gam, se.fit = TRUE, type = "response")

# Add confidence intervals
add_confidence_intervals <- function(tbl, level = 0.95) {
	z <- qnorm(1 - (1 - level) / 2)
	tbl |>
		mutate(
			lower = predicted_value - z * se,
			upper = predicted_value + z * se
		)
}


# Combine into results
BLUPs_gam_single <- newdata_gam |>
	mutate(
		predicted_value = preds_treatment_centered$fit,
		se = preds_treatment_centered$se.fit
	) |>
	add_confidence_intervals()


# Caterpillar plot
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

# Plot
BLUPs_gam_single |>
	add_confidence_intervals() |>
	plot_treatment_caterpillar()

# Residuals for model fit check
yield_data_gam_res <- yield_data_single_gam |>
	mutate(
		fitted   = predict(m2, type = "response"),
		residual = residuals(m2, type = "response")
	)

ggplot(yield_data_gam_res, aes(x = fitted, y = residual)) +
	geom_point(alpha = 0.8, color = "black") +
	geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
	labs(
		x = "Fitted values",
		y = "Residuals",
		title = "Residuals vs Fitted Values (GAM)"
	) +
	theme_minimal(base_size = 14)

# Histogram of Residuals
hist_plot <- ggplot(yield_data_gam_res, aes(x = residual)) +
	geom_histogram(binwidth = 0.3, fill = "steelblue", color = "black") +
	labs(
		title = "Histogram of Residuals (GAM)",
		x = "Residual",
		y = "Frequency"
	) +
	theme_minimal(base_size = 14)

# Residuals vs Fitted
resid_fit_plot <- ggplot(yield_data_gam_res, aes(x = fitted, y = residual)) +
	geom_point(alpha = 0.8, color = "black") +
	geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
	labs(
		x = "Fitted Values",
		y = "Residuals",
		title = "Residuals vs Fitted Values (GAM)"
	) +
	theme_minimal(base_size = 14)

# QQ Plot
qq_plot <- ggplot(yield_data_gam_res, aes(sample = residual)) +
	stat_qq(color = "darkgreen") +
	stat_qq_line(color = "black", linewidth = 0.8) +
	labs(
		title = "QQ Plot of Residuals (GAM)",
		x = "Theoretical Quantiles",
		y = "Sample Quantiles"
	) +
	theme_minimal(base_size = 14)

# Arrange in 2x2 grid (fourth empty)
(hist_plot | resid_fit_plot) /
	(qq_plot | plot_spacer())


# ASReml Single -----------------------------------------------------------
yield_data_single_asreml <- yield_data_single_gam |> 
	arrange(ROW_fc, COLUMN_fc)

library(asreml)

m2_asreml_model <- asreml(
	fixed = YIELD ~ ROW + COLUMN,
	random = ~ spl(ROW) + spl(COLUMN) + TREATMENT,
	residual = ~ ar1(ROW_fc):ar1(COLUMN_fc),
	data = yield_data_single_asreml, 
	options = options(trace=FALSE)
)

m2_asreml_model <- update(m2_asreml_model)
m2_asreml_model <- update(m2_asreml_model)

pred_asreml <- predict(m2_asreml_model, classify = "TREATMENT", sed = TRUE)$pvals
# Asreml Results
pred_asreml |> as.data.frame() |> head(10)

# GAM results
BLUPs_gam_single |> as.data.frame() |> head(10)

cor(BLUPs_gam_single$predicted_value, pred_asreml$predicted.value)

yield_data_single_asreml_res <- yield_data_single_asreml |>
	dplyr::mutate(
		fitted = fitted(m2_asreml_model),
		residual = residuals(m2_asreml_model, type = "response")
	)

hist(yield_data_single_asreml_res$residual)

ggplot(yield_data_single_asreml_res, aes(x = fitted, y = residual)) +
	geom_point(alpha = 0.8, color = "black") +
	geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
	labs(
		x = "Fitted values",
		y = "Residuals",
		title = "Residuals vs Fitted Values with LOESS Curve"
	) +
	theme_minimal(base_size = 14)

# From GAM: full model predicted means
# Use centered predictions from GAM
rank_gam_total <- BLUPs_gam_single |>
	dplyr::select(TREATMENT, predicted_value) |>
	dplyr::mutate(rank_gam_total = rank(-predicted_value))

# ASReml ranks
rank_asreml <- pred_asreml |>
	dplyr::select(TREATMENT, predicted.value) |>
	dplyr::mutate(rank_asreml = rank(-predicted.value))

# Compare ranks
rank_comparison <- rank_gam_total |>
	dplyr::inner_join(rank_asreml, by = "TREATMENT") |>
	dplyr::mutate(diff = rank_gam_total - rank_asreml)

# View largest differences
rank_comparison |> arrange(desc(abs(diff)))

cor(rank_comparison$rank_gam_total, rank_comparison$rank_asreml, method = "kendall")

ggplot(rank_comparison, aes(x = rank_asreml, y = rank_gam_total)) +
	geom_point(size = 2.5, color = "darkgreen") +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
	labs(
		title = "Rank Comparison: ASReml vs GAM (Full Fitted Values)",
		x = "ASReml Rank",
		y = "GAM Rank"
	) +
	theme_minimal(base_size = 14)

compare_top_bottom_ranks_single_env <- function(gam_df, asreml_df, n = NULL, percent = 20,
																																																reference = c("gam", "asreml")) {
	reference <- match.arg(reference)
	
	stopifnot(all(c("TREATMENT", "predicted_value") %in% names(gam_df)))
	stopifnot(all(c("TREATMENT", "predicted.value") %in% names(asreml_df)))
	
	total_genotypes <- length(unique(gam_df$TREATMENT))
	if (is.null(n)) {
		n <- ceiling(percent / 100 * total_genotypes)
		message("Selecting top and bottom ", n, " genotypes (", percent, "% of total)...")
	} else {
		message("Selecting top and bottom ", n, " genotypes ...")
	}
	
	# Rank both models
	gam_ranked <- gam_df |>
		dplyr::mutate(rank_gam = rank(-predicted_value, ties.method = "min"))
	
	asreml_ranked <- asreml_df |>
		dplyr::rename(pred_asr = predicted.value) |>
		dplyr::mutate(rank_asreml = rank(-pred_asr, ties.method = "min"))
	
	merged <- dplyr::inner_join(gam_ranked, asreml_ranked, by = "TREATMENT") |>
		dplyr::select(TREATMENT, rank_gam, rank_asreml, predicted_value, pred_asr)
	
	# Top and Bottom n from both models
	top_gam <- merged |> dplyr::arrange(rank_gam) |> dplyr::slice_head(n = n)
	top_asr <- merged |> dplyr::arrange(rank_asreml) |> dplyr::slice_head(n = n)
	bottom_gam <- merged |> dplyr::arrange(desc(rank_gam)) |> dplyr::slice_head(n = n)
	bottom_asr <- merged |> dplyr::arrange(desc(rank_asreml)) |> dplyr::slice_head(n = n)
	
	# Intersections
	top_shared <- dplyr::intersect(top_gam$TREATMENT, top_asr$TREATMENT)
	bot_shared <- dplyr::intersect(bottom_gam$TREATMENT, bottom_asr$TREATMENT)
	
	top_shared_pct <- length(top_shared) / n * 100
	bot_shared_pct <- length(bot_shared) / n * 100
	
	# Kendall tau only on shared genotypes
	top_kendall <- merged |>
		dplyr::filter(TREATMENT %in% top_shared) |>
		dplyr::summarise(
			section = "Top",
			kendall_tau = cor(rank_gam, rank_asreml, method = "kendall")
		)
	
	bottom_kendall <- merged |>
		dplyr::filter(TREATMENT %in% bot_shared) |>
		dplyr::summarise(
			section = "Bottom",
			kendall_tau = cor(rank_gam, rank_asreml, method = "kendall")
		)
	
	# Combine
	kendall_summary <- dplyr::bind_rows(top_kendall, bottom_kendall) |>
		dplyr::mutate(
			shared_pct = c(top_shared_pct, bot_shared_pct),
			section = factor(section, levels = c("Top", "Bottom"))
		) |>
		dplyr::arrange(section)
	
	# Return ranks for diagnostics
	full_subset <- dplyr::bind_rows(
		top_gam |> dplyr::mutate(section = "Top"),
		bottom_gam |> dplyr::mutate(section = "Bottom")
	)
	
	list(
		summary = kendall_summary,
		ranks = full_subset
	)
}

ggplot(yield_data_gam_res, aes(x = COLUMN, y = ROW, fill = residual)) +
	geom_tile(color = "white") +
	scale_fill_viridis(
		option = "viridis",  # blue → green → yellow
		name = "Residual",
		direction = 1
	) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	coord_fixed() +
	labs(
		title = "Spatial Heatmap of Residuals (GAM)",
		x = "Column",
		y = "Row"
	) +
	theme_minimal(base_size = 13) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5),
		legend.position = "right"
	)

ggplot(yield_data_single_asreml_res, aes(x = COLUMN, y = ROW, fill = residual)) +
	geom_tile(color = "white") +
	scale_fill_viridis(
		option = "viridis",  # blue → green → yellow
		name = "Residual",
		direction = 1
	) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	coord_fixed() +
	labs(
		title = "Spatial Heatmap of Residuals (ASReml)",
		x = "Column",
		y = "Row"
	) +
	theme_minimal(base_size = 13) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5),
		legend.position = "right"
	)


results <- compare_top_bottom_ranks_single_env(
	BLUPs_gam_single, 
	pred_asreml,
	n = 20,
	reference = "asreml"
)

results

rankings <- data.frame(
	Top = c(10, 20, 30, 40, 50),
	kendall_tau = c(1.0, 1.0, 0.991, 0.987, 0.982),
	shared_pct = c(100, 100, 100, 97.5, 100)
) |> kable()



library(dplyr)
library(tidyr)
library(ggplot2)

# Prepare long-format data
se_long <- BLUPs_gam_single |>
	select(TREATMENT, se_gam = se) |>
	left_join(
		pred_asreml |> select(TREATMENT, se_asreml = std.error),
		by = "TREATMENT"
	) |>
	pivot_longer(
		cols = c(se_gam, se_asreml),
		names_to = "model",
		values_to = "se"
	) |>
	mutate(
		model = recode(model, se_gam = "GAM", se_asreml = "ASReml"),
		TREATMENT = factor(TREATMENT, levels = unique(TREATMENT))  # preserve order
	)

# Plot
ggplot(se_long, aes(x = TREATMENT, y = se, group = model, color = model)) +
	geom_line(linewidth = 0.8) +
	labs(
		title = "Standard Errors of Treatment Means",
		x = "Treatment",
		y = "Standard Error",
		color = "Model"
	) +
	theme_classic(base_size = 14) +
	theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
		plot.title = element_text(face = "bold")
	)