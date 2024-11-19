import os

from flask import Flask, request, jsonify, send_file
from flask_cors import cross_origin
import polars as pl
import io

app = Flask(__name__)

DATA_PATH = os.environ.get("BATCH_QUERY_DATA_PATH", "./batch_query_data_parquet")
dataset = pl.read_parquet(f"{DATA_PATH}/*.parquet")


def print_phenotype(phenotype):
    return (
        f'id: {phenotype["id"]}, name: {phenotype["name"]}'
        if phenotype is not None
        else ""
    )


def print_phenotype_list(phenotype_list):
    return phenotype_list.map_elements(print_phenotype, return_dtype=pl.Utf8)


def flatten_nested_columns(input_df):
    result_df = input_df.with_columns(
        displayPhenotype=pl.when(pl.col("displayPhenotype").is_not_null()).then(
            pl.col("displayPhenotype").map_elements(
                print_phenotype, return_dtype=pl.Utf8
            )
        ),
        significantPhenotype=pl.when(pl.col("significantPhenotype").is_not_null()).then(
            pl.col("significantPhenotype").map_elements(
                print_phenotype, return_dtype=pl.Utf8
            )
        ),
        phenotypeSexes=pl.when(pl.col("phenotypeSexes").is_not_null()).then(
            pl.col("phenotypeSexes").cast(pl.List(pl.Utf8)).list.join(", ")
        ),
        stringifiedList_intermediatePhenotypes=pl.when(
            pl.col("intermediatePhenotypes").is_not_null()
        ).then(
            pl.col("intermediatePhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ),
        stringifiedList_potentialPhenotypes=pl.when(
            pl.col("potentialPhenotypes").is_not_null()
        ).then(
            pl.col("potentialPhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ),
        stringifiedList_topLevelPhenotypes=pl.when(
            pl.col("topLevelPhenotypes").is_not_null()
        ).then(
            pl.col("topLevelPhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ),
        stringifiedList_humanPhenotypes=pl.when(
            pl.col("humanPhenotypes").is_not_null()
        ).then(
            pl.col("humanPhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ),
    )

    result_df = result_df.drop("intermediatePhenotypes")
    result_df = result_df.drop("potentialPhenotypes")
    result_df = result_df.drop("topLevelPhenotypes")
    result_df = result_df.drop("humanPhenotypes")

    result_df = result_df.rename(
        {
            "stringifiedList_intermediatePhenotypes": "intermediatePhenotypes",
            "stringifiedList_potentialPhenotypes": "potentialPhenotypes",
            "stringifiedList_topLevelPhenotypes": "topLevelPhenotypes",
            "stringifiedList_humanPhenotypes": "humanPhenotypes",
        }
    )

    return result_df


@app.route("/mi/impc/batch-query", methods=["POST"])
@cross_origin()
def query_data():
    # Parse request headers and data
    response_format = request.headers.get("Accept", "application/json").lower()
    mgi_ids = []

    if "file" in request.files:
        file = request.files["file"]
        mgi_ids = file.read().decode("utf-8").strip().split("\n")
    elif request.json and "mgi_ids" in request.json:
        mgi_ids = request.json["mgi_ids"]
    else:
        return jsonify({"error": "No MGI accession IDs provided"}), 400

    # Query the dataset
    filtered_data = dataset.filter(pl.col("mgiGeneAccessionId").is_in(mgi_ids))

    if response_format == "application/json":
        return jsonify(filtered_data.to_dicts())
    elif response_format == "xlsx":
        return dataframe_to_xlsx(filtered_data)
    elif response_format == "tsv":
        return dataframe_to_tsv(filtered_data)
    else:
        return jsonify({"error": "Unsupported response format"}), 400


def dataframe_to_xlsx(df):
    output = io.BytesIO()
    new_df = flatten_nested_columns(df)
    new_df.write_excel(output)
    output.seek(0)
    return send_file(
        output,
        mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        as_attachment=True,
        download_name="query_results.xlsx",
    )


def dataframe_to_tsv(df):
    output = io.BytesIO()
    new_df = flatten_nested_columns(df)
    new_df.write_csv(output, separator="\t")
    output.seek(0)
    return send_file(
        output,
        mimetype="text/tab-separated-values",
        as_attachment=True,
        download_name="query_results.tsv",
    )


@app.route("/mi/impc/batch-query/health-check", methods=["GET"])
def health_check():
    return jsonify({"status": "ok"}), 200


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
