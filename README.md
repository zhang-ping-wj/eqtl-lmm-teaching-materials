# Teaching note: a simulated eQTL analysis with repeated samples

This note is based on the files in `COI_summer.teaching/`:

- `Simulated_test_data_geneA_eQTL.csv`
- `Simulated test dataset_geneA_eQTL.R`

The workshop document describes a practical on eQTLs, genotype coding, linear mixed-effects models, and genotype-by-time analysis. The CSV and R script are simulated teaching data, not a genome-wide eQTL mapping.

## Data structure summary

The dataset is in **long format**: each row is one sample from one individual at one time point. There are 510 sample rows and 300 unique individuals (`subject_id`), so the rows are not all independent observations.

| Feature | Observed in the files |
|---|---:|
| Individuals | 300 |
| Samples/rows | 510 |
| Baseline samples | 300 |
| 6-hour post-infection samples | 62 |
| 24-hour post-infection samples | 148 |
| Individuals with baseline only | 152 |
| Individuals with baseline + 24 hours | 86 |
| Individuals with all three time points | 62 |
| Genotype groups (0, 1, 2 effect-allele copies) | 270, 203, 37 samples |
| Sex | F: 246; M: 264 |
| Batch | Batch1: 211; Batch2: 168; Batch3: 131 |

Available columns are `subject_id`, `sample_id`, `genotype_effect_allele_dosage`, `age_years`, `sex`, `batch`, `timepoint`, `hours_post_infection`, `infection_context`, and `geneA_expression_log2`.

The outcome appears to be log2 Gene A expression. Predictors include genotype dosage, age, sex, batch, and infection time/context. Genotype is constant within an individual, while expression and time vary by sample. This is a longitudinal/repeated-measures dataset with incomplete follow-up: all 300 individuals have baseline, but only 62 have all three time points.

## Why use an LMM here?

A simple linear model treats every row as if it came from a different, unrelated individual. That is not true here: some people contribute measurements at several time points. Samples from one person tend to be more alike because they share biology, genetics, and other unmeasured characteristics.

The script uses models such as:

```r
lmer(expression ~ genotype + age + sex + batch + (1 | subject_id), data = df)
```

The fixed effects estimate average relationships across the sample. The random intercept `(1 | subject_id)` gives each individual their own baseline level, capturing stable person-to-person differences and the correlation among that person's repeated samples. Without it, standard errors could be too small because repeated observations would be counted as if they were entirely new people.

The random effect is not the genotype effect. Genotype is a predictor whose average effect is estimated directly. The random intercept captures unexplained individual-specific starting levels. A more elaborate analysis could consider a random slope for time, but that is not fitted here.

The script first adjusts expression for age, sex, batch, and the subject random intercept, then defines `corrected_exp` as residual plus the overall mean. This is intended as a teaching simplification; the subsequent interaction model is fitted to this corrected outcome rather than writing all covariates in one formula.

## How to read the volcano plot

Important file-specific caveat: **the supplied R script does not create a conventional volcano plot**. It creates a genotype violin plot with beta and p-value annotation, followed by a genotype-by-time interaction plot. The workshop document mentions a volcano plot for a separate ATAC-seq/ChIP-seq practical, but no volcano-plot data or code is present in this folder.

If a conventional eQTL volcano plot is introduced as an extension:

- The x-axis is usually the estimated genetic effect, such as a regression coefficient or log2 fold change. Right means higher expression with increasing effect-allele dosage; left means lower expression.
- The y-axis is usually `-log10(p-value)`. Higher points have smaller p-values and stronger statistical evidence against no association.
- Each point represents one tested variant-gene association.

This dataset has one named outcome, Gene A expression, and one genotype dosage, so it supports a single-association demonstration rather than a full cloud of genome-wide points. Students should separate effect size (how large and in which direction) from statistical evidence (how compatible the data are with zero effect). A point high near x = 0 can be statistically convincing but biologically small; a point far from zero but low can be imprecisely estimated. In a real multi-test eQTL scan, multiple-testing correction would be essential. The supplied script does not perform it.

For the actual genotype plot, ask whether expression distributions differ across dosage groups 0, 1, and 2, and whether the beta and p-value support that pattern. Violin shapes show distributions, jittered points show samples, and the boxplot shows median and spread.

### Generated plot: genotype violin plot

![Gene A expression by genotype](geneA_violin_plot.png)

*Figure 1. Simulated Gene A expression after adjustment for age, sex, batch, and subject-level differences, shown across 0, 1, and 2 copies of the effect allele. The displayed beta is the fitted average change per additional effect-allele copy.*

## How to read the interaction plot

The interaction model is:

```r
lmer(corrected_exp ~ genotype_effect_allele_dosage * timepoint +
       (1 | subject_id), data = df)
```

The `*` includes genotype, time point, and genotype-by-time point terms. The interaction asks whether the genotype slope is different at different time points.

In the plot, the x-axis is genotype dosage, the y-axis is corrected expression, and each coloured line is a time point. The overall interaction test compares a model with genotype and time point only against one that also allows genotype effects to vary by time.

- Roughly parallel lines suggest a similar genotype effect over time, even if expression changes with time.
- Different slopes suggest a time-dependent or context-dependent eQTL.
- Converging or crossing lines suggest that infection changes the difference between genotype groups, or even their ordering.
- A vertical shift with parallel lines suggests a time main effect, not necessarily an interaction.

Biologically, a time-dependent slope could mean infection changes the regulatory environment around Gene A, for example signalling, chromatin state, transcription-factor activity, or cell composition. The plot alone does not establish mechanism or causality.

Because the data are unbalanced, some genotype-by-time combinations are small, especially genotype 2 and 6 hours. Students should inspect points and uncertainty intervals, not only fitted lines. Time is treated as a factor here, so the model estimates separate categorical differences rather than assuming a straight-line change from 0 to 6 to 24 hours.

### Generated plot: genotype-by-time interaction

![Gene A genotype by infection time interaction](geneA_interaction_plot.png)

*Figure 2. Fitted genotype-expression relationships at baseline, 6 hours, and 24 hours after infection. The non-parallel lines illustrate a time-dependent eQTL effect in these simulated data; shaded bands show model-based uncertainty.*

## File-specific assumptions and limitations

- The data are labelled simulated, so biological conclusions are illustrative only.
- The R script reads from `Downloads/Simulated_test_data_geneA_eQTL.csv`, whereas the supplied CSV is in the teaching folder; the path may need changing.
- The first output is a violin plot, not a volcano plot. A conventional volcano plot requires results for many tested associations, which are not included.
- The analysis uses a subject random intercept but no random time slope, explicit residual correlation structure, or multiple-testing correction.
- Follow-up is incomplete: all 300 individuals have baseline, but not all have both post-infection samples.

## Suggested figure captions

**Genotype violin plot:** “Simulated Gene A expression, adjusted for age, sex, batch, and subject-level baseline differences, shown across 0, 1, and 2 copies of the effect allele. Points are samples; violins show distributions; the annotation reports the fitted genotype effect and p-value.”

**Interaction plot:** “Simulated genotype-by-infection-time analysis for Gene A. Coloured fitted lines show the estimated relationship between effect-allele dosage and corrected expression at baseline, 6 hours, and 24 hours after infection. Non-parallel slopes indicate that the eQTL effect may depend on infection time.”
