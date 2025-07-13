# import basic packages
import pandas as pd
from snakemake.utils import validate


# read sample sheet
dnaseq = (
    pd.read_csv(config["dnaseq_samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

scrnaseq = (
    pd.read_csv(config["scrnaseq_samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("id", drop=False)
    .sort_index()
)

# validate sample sheet and config file
validate(dnaseq, schema="../../config/schemas/dnaseq.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")
