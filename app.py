from flask import Flask, request, jsonify, send_file
import polars as pl
import io

app = Flask(__name__)

# Load the dataset (assuming it's already stored in the same directory)
dataset = pl.read_parquet('batch_query_service_data_parquet/*.parquet')


@app.route('/query', methods=['POST'])
def query_data():
    # Parse request headers and data
    response_format = request.headers.get('Response-Format', 'json').lower()
    mgi_ids = []

    if 'file' in request.files:
        file = request.files['file']
        mgi_ids = file.read().decode('utf-8').strip().split('\n')
    elif request.json and 'mgi_ids' in request.json:
        mgi_ids = request.json['mgi_ids']
    else:
        return jsonify({"error": "No MGI accession IDs provided"}), 400

    # Query the dataset
    filtered_data = dataset.filter(pl.col('mgiGeneAccessionId').is_in(mgi_ids))

    if response_format == 'json':
        return jsonify(filtered_data.to_dicts())
    elif response_format == 'xlsx':
        return dataframe_to_xlsx(filtered_data)
    elif response_format == 'tsv':
        return dataframe_to_tsv(filtered_data)
    else:
        return jsonify({"error": "Unsupported response format"}), 400


def dataframe_to_xlsx(df):
    output = io.BytesIO()
    df.write_excel(output)
    output.seek(0)
    return send_file(output, mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     as_attachment=True, download_name='query_results.xlsx')


def dataframe_to_tsv(df):
    output = io.BytesIO()
    df.write_csv(output, separator='\t')
    output.seek(0)
    return send_file(output, mimetype='text/tab-separated-values', as_attachment=True,
                     download_name='query_results.tsv')


if __name__ == '__main__':
    app.run(debug=True)
