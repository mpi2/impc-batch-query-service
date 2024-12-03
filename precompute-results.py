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


def get_filtered_list(df, alleleSymbol, dfColumn, significant_data=True):
    filteredDf = df.filter((pl.col("alleleSymbol") == alleleSymbol) & (
        pl.col("significant") == significant_data))
    list_with_data = filteredDf[dfColumn].unique().to_list()
    return list(filter(None, list_with_data))


def get_filtered_significant_systems(df, alleleSymbol, significant_data=True):
    filteredDf = df.filter((pl.col("alleleSymbol") ==
                           alleleSymbol) & (pl.col("significant") == significant_data))
    filteredDf = filteredDf.with_columns(
        systemNames=pl.concat_list("topLevelPhenotypeNames")
    )
    return list(dict.fromkeys(sum(filter(None, filteredDf["systemNames"].to_list()), [])))


def get_phenotype_names(phenotype_list):
    return phenotype_list.map_elements(lambda p: p["name"], return_dtype=pl.Utf8)


def ensemble_allele_result(allele, df):
    significantLifeStages = get_filtered_list(df, allele, "lifeStageName")
    significantPhenotypes = get_filtered_list(
        df, allele, "significantPhenotypeName")
    notSignificantPhenotypes = get_filtered_list(
        df, allele, "significantPhenotypeName", False)
    significantSystems = get_filtered_significant_systems(df, allele)
    notSignificantSystems = get_filtered_significant_systems(df, allele, False)

    significantLifeStages.sort()
    significantPhenotypes.sort()
    significantSystems.sort()
    notSignificantPhenotypes.sort()
    notSignificantSystems.sort()

    alleleData = {
        "allele": allele,
        "significantLifeStages": significantLifeStages,
        "significantPhenotypes": significantPhenotypes,
        "notSignificantPhenotypes": notSignificantPhenotypes,
        "significantSystems": significantSystems,
        "notSignificantSystems": notSignificantSystems
    }

    return alleleData


def get_all_significant_phenotypes(df):
    filteredDf = df.filter(pl.col("significant") == True)
    list_with_data = filteredDf["significantPhenotypeName"].unique().to_list()
    return list(filter(None, list_with_data))


def get_all_significant_systems(df):
    filteredDf = df.filter(pl.col("significant") == True)
    filteredDf = filteredDf.with_columns(
        systemNames=pl.concat_list("topLevelPhenotypeNames")
    )
    return list(dict.fromkeys(sum(filter(None, filteredDf["systemNames"].to_list()), [])))


def ensemble_gene_result(df):
    human_gene_symbols = get_value_from_df(df, "humanGeneSymbol")
    human_gene_ids = get_value_from_df(df, "hgncGeneAccessionId")
    gene_id = get_value_from_df(df, "id")
    alleles = get_value_from_df(df, "alleleSymbol", unwrap_value=False)
    mouse_gene_symbol = df["alleleSymbol"].to_list()[0].split("<")[0]
    all_significant_systems = get_all_significant_systems(df)
    all_significant_phenotypes = get_all_significant_phenotypes(df)

    alelles_data = list(map(functools.partial(
        ensemble_allele_result, df=df), alleles))

    geneData = {
        "mouseGeneSymbol": mouse_gene_symbol,
        "allSignificantSystems": all_significant_systems,
        "allSignificantPhenotypes": all_significant_phenotypes,
        "humanGeneSymbols": human_gene_symbols,
        "humanGeneIds": human_gene_ids,
        "geneId": gene_id,
        "alleles": alelles_data,
    }
    return geneData


def process_gene(full_dataset, mgi_id):
    gene_data = full_dataset.filter(pl.col("mgiGeneAccessionId").eq(mgi_id))
    df = gene_data.drop(["statisticalResultId", "statisticalResultId", "alleleName",
                        "dataType", "effectSize", "femaleMutantCount", "maleMutantCount", "metadataGroup", "pValue", "parameterName",
                         "parameterStableId", "phenotypeSexes", "phenotypingCentre", "pipelineStableId", "procedureMinAnimals",
                         "procedureMinFemales", "procedureMinMales", "procedureName", "procedureStableId", "projectName",
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
        print(
            f"\rprocessing GENE {mgi_id}, {index} of {ids_count}", end="", flush=True)
        full_results[mgi_id] = process_gene(full_dataset, mgi_id)

    with open("preprocessed-results.json", "w") as outfile:
        print("\nDONE")
        json.dump(full_results, outfile)


if __name__ == "__main__":
    main()
