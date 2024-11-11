import os
import polars as pl
import functools
import json

DATA_PATH = os.environ.get("BATCH_QUERY_DATA_PATH", ".")


def get_value_from_df(df, dfKey, unwrap_value=True):
    val = df[dfKey].unique().to_list()
    if len(val) == 1 and unwrap_value:
        val = val[0]
    return val


def get_filtered_list(df, alleleSymbol, dfColumn):
    filteredDf = df.filter(pl.col("alleleSymbol") == alleleSymbol)
    list_with_data = filteredDf[dfColumn].unique().to_list()
    return list(filter(None, list_with_data))


def get_phenotype_names(phenotype_list):
    return phenotype_list.map_elements(lambda p: p["name"], return_dtype=pl.Utf8)


def ensemble_allele_result(allele, df):
    df = df.with_columns(
        mergedTopLevelPhenotypeNames=pl.concat_list("topLevelPhenotypeNames")
    )
    significantLifeStages = get_filtered_list(df, allele, "lifeStageName")
    significantPhenotypes = get_filtered_list(
        df, allele, "significantPhenotypeName")
    significantSystems = list(dict.fromkeys(sum(
        filter(None, df["mergedTopLevelPhenotypeNames"].to_list()), [])))

    significantLifeStages.sort()
    significantPhenotypes.sort()
    significantSystems.sort()

    alleleData = {
        "allele": allele,
        "significantLifeStages": significantLifeStages,
        "significantPhenotypes": significantPhenotypes,
        "significantSystems": significantSystems
    }

    return alleleData


def ensemble_gene_result(df):
    humanGeneSymbols = get_value_from_df(df, "humanGeneSymbol")
    humanGeneIds = get_value_from_df(df, "hgncGeneAccessionId")
    geneId = get_value_from_df(df, "id")
    alleles = get_value_from_df(df, "alleleSymbol", unwrap_value=False)
    alellesData = list(map(functools.partial(
        ensemble_allele_result, df=df), alleles))

    geneData = {
        "humanGeneSymbols": humanGeneSymbols,
        "humanGeneIds": humanGeneIds,
        "geneId": geneId,
        "alleles": alellesData,
    }
    return geneData


def process_gene(full_dataset, mgi_id):
    gene_data = full_dataset.filter(pl.col("mgiGeneAccessionId").eq(mgi_id))
    df = gene_data.drop(["statisticalResultId", "statisticalResultId", "alleleName",
                        "dataType", "effectSize", "femaleMutantCount", "maleMutantCount", "metadataGroup", "pValue", "parameterName",
                         "parameterStableId", "phenotypeSexes", "phenotypingCentre", "pipelineStableId", "procedureMinAnimals",
                         "procedureMinFemales", "procedureMinMales", "procedureName", "procedureStableId", "projectName", "significant",
                         "statisticalMethod", "status", "zygosity", "intermediatePhenotypes", "potentialPhenotypes", "humanPhenotypes",
                         "alleleAccessionId"
                         ])

    df = df.with_columns(
        displayPhenotypeName=pl.when(
            pl.col("displayPhenotype").is_not_null()
        )
        .then(
            pl.col("displayPhenotype").map_elements(
                lambda phenotype: phenotype["name"], return_dtype=pl.Utf8)
        ),
        displayPhenotypeId=pl.when(
            pl.col("displayPhenotype").is_not_null()
        )
        .then(
            pl.col("displayPhenotype").map_elements(
                lambda phenotype: phenotype["id"], return_dtype=pl.Utf8)
        ),
        significantPhenotypeName=pl.when(
            pl.col("significantPhenotype").is_not_null()
        )
        .then(
            pl.col("significantPhenotype").map_elements(
                lambda phenotype: phenotype["name"], return_dtype=pl.Utf8)
        ),
        significantPhenotypeId=pl.when(
            pl.col("significantPhenotype").is_not_null()
        )
        .then(
            pl.col("significantPhenotype").map_elements(
                lambda phenotype: phenotype["id"], return_dtype=pl.Utf8)
        ),
        topLevelPhenotypeNames=pl.when(
            pl.col("topLevelPhenotypes").is_not_null()
        )
        .then(
            pl.col("topLevelPhenotypes")
            .map_elements(get_phenotype_names, return_dtype=pl.List(pl.Utf8))
        )
    )
    df = df.drop(
        ["displayPhenotype", "topLevelPhenotypes", "significantPhenotype"])
    df = df.unique()
    df = df.with_row_index()

    return ensemble_gene_result(df)


def create_dataset():
    return pl.read_parquet(f"{DATA_PATH}/*.parquet")


def main():
    full_dataset = create_dataset()
    mgi_ids = full_dataset["mgiGeneAccessionId"].unique().to_list()
    mgi_ids = list(filter(None, mgi_ids))
    ids_count = len(mgi_ids)
    print(f"total of {ids_count} genes")
    full_results = {}

    for index, mgi_id in enumerate(mgi_ids):
        print(f"processing GENE {mgi_id}, {index} of {ids_count}")
        full_results[mgi_id] = process_gene(full_dataset, mgi_id)

    with open("full-dataset-results.json", "w") as outfile:
        json.dump(full_results, outfile)


if __name__ == "__main__":
    main()
