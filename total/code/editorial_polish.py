"""One-time, reversible editorial pass; uses existing results only."""
from pathlib import Path
import shutil, json
from datetime import datetime, timezone

root = Path(__file__).resolve().parents[1]
run = root / 'reproducibility' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_editorial')
run.mkdir(parents=True)
names = ['build_manuscript.py', 'build_supplementary.py', 'fig6_mechanism.py', 'fig2_development.py']
for name in names:
    shutil.copy2(root / 'code' / name, run / name)
for name in ['ZEB1_TALL_manuscript.docx', 'ZEB1_TALL_Supplementary.docx']:
    shutil.copy2(root / 'manuscript' / name, run / name)

p = root / 'code/build_manuscript.py'
s = p.read_text(encoding='utf-8')
s = s.replace('"ZEB1 expression in T-cell acute lymphoblastic leukemia varies with "\n              "developmental and molecular subtype context: "\n              "a systematic multi-cohort, multi-omic reappraisal"', '"Context-dependent ZEB1 expression and treatment-response associations "\n              "in T-cell acute lymphoblastic leukemia"')
s = s.replace('Word count (main text): ~5,400  |  Main figures: 7  |  Supplementary figures: 15  |  Cohorts: >15 public datasets', 'Main figures: 7  |  Supplementary figures: 15')
start = s.index('body(\n', s.index('H("Abstract"'))
end = s.index('\n# ===========================================================================', start)
abstract = ('Developmental identity may shape how ZEB1 relates to leukemic cell state and treatment response in T-cell acute lymphoblastic leukemia (T-ALL). '
    'We integrated public bulk and single-cell transcriptomes, chromatin profiles, perturbations, drug screens and clinical endpoints to define these relationships across biological contexts. '
    'Across three cohorts, ZEB1 expression was higher in T-ALL than B-ALL (pooled Hedges g = 1.44). Normal thymocyte references showed developmental-stage variation, while TLX-high T-ALL showed directionally higher ZEB1 across four cohorts (pooled g = 0.61; the small-study-adjusted interval included zero). '
    'In a 108-patient NOPHO cohort, ZEB1 correlated positively with BCL11B and IL7R and inversely with LMO2. Perturbation responses depended on context: TLX1 knockdown increased ZEB1, and Lmo2-knockout single-cell data showed higher Zeb1 at day 3 under two priming conditions. Patient single-cell associations varied with cohort and malignant-cell selection. '
    'Higher ZEB1 expression correlated with higher methotrexate IC50 in GDSC2 (nine T-ALL models; Spearman rho = 0.90) and GDSC1 (ten models; rho = 0.67); eight shared models limit independence. CRISPR effects were heterogeneous, with four of eight T-ALL models crossing the conventional gene-effect reference of -0.5. '
    'In patient data, none of 17 analyzable ex-vivo drug correlations met FDR < 0.05; methotrexate was not assayed. Adding diagnosis ZEB1 to mid-induction residual disease did not improve internal prediction of later positivity in an exploratory 66-patient subset. '
    'These findings connect ZEB1 to developmental and molecular context and identify a treatment-specific association with methotrexate response in model systems. Clinical predictive value remains unresolved, motivating context-defined interpretation of ZEB1 biology.')
s = s[:start] + 'body(\n    ' + repr(abstract) + '\n)\n' + s[end:]
s = s.replace('Epigenomic and perturbation data give correlative, not causal, support', 'Chromatin context and perturbation responses vary across developmental models')
s = s.replace('Mechanistic layers are correlative. (A)', 'Chromatin context and transcriptional responses to perturbation. (A)')
s = s.replace('Available pharmacogenomic data do not support a general ZEB1 therapeutic biomarker', 'ZEB1 associations with drug response and genetic dependency vary by context')
s = s.replace('"Available pharmacogenomic data do not support a general ZEB1 therapeutic "\n       "biomarker. (A)', '"Drug response, genetic dependency and residual disease. (A)')
s = s.replace('There is no independent "\n       "clinical validation was available.', 'Independent "\n       "clinical validation was not available.')
start = s.index('body(\n', s.index('H("Discussion"'))
end = s.index('body(\n', start + 6)
opening = ('ZEB1 expression in T-ALL is embedded in developmental and molecular context. Its higher expression relative to B-ALL, variation across normal thymocyte stages and positive TLX-high direction across cohorts provide a framework for interpreting heterogeneous patient and model-system findings. The positive BCL11B and IL7R associations and inverse LMO2 association in NOPHO, together with context-specific perturbation responses, support a biological relationship worth resolving at the level of defined cell states. Bulk associations and selected model systems do not by themselves identify a universal regulatory mechanism; patient single-cell analyses further show that malignant-cell selection can alter the apparent relationship.')
clinical = ('The clinically relevant question is whether this context helps distinguish treatment response. The methotrexate association is a focused lead: higher ZEB1 tracked higher IC50 in both GDSC screens, although eight overlapping models prevent an independent-replication claim. This direction differs from a general ZEB1-low chemoresistance model and suggests that drug identity, developmental state and model background should be considered together. The patient ex-vivo panel did not contain methotrexate, so its results cannot directly adjudicate this signal. Separately, diagnosis ZEB1 did not improve internal prediction of persistent MRD beyond mid-induction MRD in the available subset. Drug sensitivity in cultured models and incremental patient-level prediction address different questions; the present data connect the former to a specific treatment while leaving clinical translation unresolved.')
s = s[:start] + 'body(\n    ' + repr(opening) + '\n)\nbody(\n    ' + repr(clinical) + '\n)\n' + s[end:]
s = s.replace('These conclusions refine, rather than simply negate, prior work. The TLX-high ZEB1 ', 'The TLX-high ZEB1 ')
s = s.replace('    p.paragraph_format.space_before = Pt(before)\n', '    p.paragraph_format.keep_with_next = True\n    p.paragraph_format.space_before = Pt(before)\n')
s = s.replace('    fmt = p.paragraph_format\n', '    fmt = p.paragraph_format\n    fmt.widow_control = True\n')
# Keep figure image and the first caption lines together; allow long legends to flow.
s = s.replace('    p.paragraph_format.space_before = Pt(6)\n', '    p.paragraph_format.keep_with_next = True\n    p.paragraph_format.space_before = Pt(6)\n')
p.write_text(s, encoding='utf-8')

p = root / 'code/build_supplementary.py'
s = p.read_text(encoding='utf-8')
s = s.replace('"ZEB1 expression in T-cell acute lymphoblastic leukemia varies with "\n              "developmental and molecular subtype context"', '"Context-dependent ZEB1 expression and treatment-response associations "\n              "in T-cell acute lymphoblastic leukemia"')
s = s.replace('    p.paragraph_format.space_before = Pt(before)\n', '    p.paragraph_format.keep_with_next = True\n    p.paragraph_format.space_before = Pt(before)\n')
s = s.replace('Hallmark/Reactome/GO are null or artefactual', 'Pathway enrichment and sensitivity to the broad expression axis')
s = s.replace('One-line claim', 'Evidence summary')
p.write_text(s, encoding='utf-8')

for name, replacements in {
    'fig6_mechanism.py': {'Promoter CpGs stay hypomethylated': 'Promoter CpG methylation', 'Shared promoter-region ChIP signal': 'Shared-region ChIP signal'},
    'fig2_development.py': {'Independent thymus atlas (Park/HTA): ZEB1 is not a late-SP gene  |  each point = one donor, n_cells≥20': 'ZEB1 across thymocyte stages | independent Park/HTA atlas', 'replicates FACS scRNA (panel B): DN/DP > SP': 'One point per donor-stage group; at least 20 cells'},
}.items():
    p = root / 'code' / name
    s = p.read_text(encoding='utf-8')
    for old, new in replacements.items():
        assert old in s, old
        s = s.replace(old, new)
    p.write_text(s, encoding='utf-8')
(run / 'scope.json').write_text(json.dumps({'scope': 'Editorial and layout pass using existing analyses', 'backup': str(run), 'abstract_words': len(abstract.split()), 'changed_code': names}, indent=2), encoding='utf-8')
print(run)
print('Abstract words:', len(abstract.split()))
