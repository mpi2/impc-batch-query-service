from flask import Flask, request, jsonify, send_file
import polars as pl
import io

app = Flask(__name__)

# Load the dataset (assuming it's already stored in the same directory)
dataset = pl.read_parquet("batch_query_service_data_parquet/*.parquet")


def print_phenotype(phenotype):
    return f'id: {phenotype["id"]}, name: {phenotype["name"]}'

 
def print_phenotype_list(phenotype_list):
    return phenotype_list.map_elements(print_phenotype, return_dtype=pl.Utf8)


def flatten_nested_columns(input_df):
    result_df = input_df
    result_df = result_df.with_columns(
        pl.col("displayPhenotype").map_elements(print_phenotype, return_dtype=pl.Utf8)
    )
    result_df = result_df.with_columns(
        pl.col("significantPhenotype").map_elements(
            print_phenotype, return_dtype=pl.Utf8
        )
    )
    result_df = result_df.with_columns(
        pl.col("phenotypeSexes")
        .cast(pl.List(pl.Utf8))
        .list.join(", ")
        .alias("stringifiedList_phenotypeSexes")
    )

    result_df = result_df.with_columns(
        (
            pl.col("intermediatePhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ).alias("stringifiedList_intermediatePhenotypes")
    )
    result_df = result_df.with_columns(
        (
            pl.col("potentialPhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ).alias("stringifiedList_potentialPhenotypes")
    )
    result_df = result_df.with_columns(
        (
            pl.col("topLevelPhenotypes")
            .map_elements(print_phenotype_list, return_dtype=pl.List(pl.Utf8))
            .list.join(" | ")
        ).alias("stringifiedList_topLevelPhenotypes")
    )

    result_df = result_df.drop("phenotypeSexes")
    result_df = result_df.drop("intermediatePhenotypes")
    result_df = result_df.drop("potentialPhenotypes")
    result_df = result_df.drop("topLevelPhenotypes")

    result_df = result_df.rename(
        {
            "stringifiedList_phenotypeSexes": "phenotypeSexes",
            "stringifiedList_intermediatePhenotypes": "intermediatePhenotypes",
            "stringifiedList_potentialPhenotypes": "potentialPhenotypes",
            "stringifiedList_topLevelPhenotypes": "topLevelPhenotypes",
        }
    )

    return result_df


@app.route("/mi/impc/batch-query", methods=["POST"])
def query_data():
    # Parse request headers and data
    response_format = request.headers.get("Response-Format", "json").lower()
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

    if response_format == "json":
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
    app.run(debug=True)
