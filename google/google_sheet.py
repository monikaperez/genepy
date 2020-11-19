from __future__ import print_function

import os
from typing import List, Optional, Tuple
from collections import defaultdict
from google.oauth2 import service_account
from googleapiclient.discovery import build
import gspread

SCOPES = [
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive.metadata.readonly",
]


def colnum_string(n):
    """
    https://stackoverflow.com/questions/23861680/convert-spreadsheet-number-to-column-letter
    Assumes n is 1-indexed
    """
    string = ""
    while True:
        n, remainder = divmod(n - 1, 26)
        string = chr(ord("A") + remainder) + string
        if n == 0:
            break
    return string


def get_cell_from_index(row_idx: int, col_idx: int):
    return colnum_string(col_idx + 1) + str(row_idx + 1)


def get_credentials():
    # return gspread.login('jkalfon@broadinstitute.org', 'JEREM1//kal')
    #gc = gspread.oauth()
    # Fetch credentials from storage
    credentials = STORAGE.get()
    # If the credentials doesn't exist in the storage location then run the flow
    if credentials is None or credentials.invalid:
        flow = flow_from_clientsecrets(CLIENT_SECRET, scope=SCOPE)
        http = httplib2.Http()
        credentials = run_flow(flow, STORAGE, http=http)
    return credentials
    credentials = authorize_credentials()


import gspread
from oauth2client.service_account import ServiceAccountCredentials


class GSheet:
    def __init__(self, creds: str, document_id: str, sheet_id: str):
        #os.system('mkdir ~/.config/')
        os.system('mkdir ~/.config/gspread/')
        os.system('mv ' + creds + " ~/.config/gspread/credentials.json")
        self.credientials = get_credentials()
        self.document_id = document_id
        self.sheet_id = sheet_id
        self.size_cache = {}
        self.cache = {}
        self.sheet_name_cache = {}
        self.request_counts = defaultdict(lambda: {"read": 0, "write": 0})
        self.client = gspread.authorize(self.credentials)

    def upload_pandas(self, df):
        df.to_csv('/tmp/googlesheetfunctiontopanda.csv')
        self.upload_csv('/tmp/googlesheetfunctiontopanda.csv')

    def upload_csv(self, file):
        spreadsheet = self.client.open(self.document_id)

        with open(file, 'r') as file_obj:
            content = file_obj.read()
            client.import_csv(spreadsheet.id, data=content)

    def get_sheet_url(self):
        return f"https://docs.google.com/spreadsheets/d/{self.document_id}/edit#gid={self.sheet_id}"

    def get_last_modified_date(self):
        response = (
            self.client.files()
            .get(fileId=self.document_id, fields="modifiedTime")
            .execute()
        )
        self.request_counts[self.document_id]["read"] += 1
        return response["modifiedTime"]

    def get_annotation_cell(self, annotation) -> Optional[str]:
        sheet = self.read_sheet()
        for row_idx, row in enumerate(sheet):
            for col_idx, val in enumerate(row):
                if val == annotation:
                    return get_cell_from_index(row_idx, col_idx)
        return None

    def _get_sheet_by_id(self, sheets, id):
        matching_sheets = [
            sheet for sheet in sheets if str(sheet["properties"]["sheetId"]) == id
        ]
        assert len(matching_sheets) == 1, f"sheets: {matching_sheets}"
        return matching_sheets[0]

    # Gets the upper bound of size (no ragged right/down)
    def get_size(self):
        key = (self.document_id, self.sheet_id)
        if key in self.size_cache:
            return self.size_cache[key]

        s = self.client.spreadsheets().get(spreadsheetId=self.document_id).execute()
        self.request_counts[self.document_id]["read"] += 1
        sheet = self._get_sheet_by_id(s["sheets"], self.sheet_id)
        gridProperties = sheet["properties"]["gridProperties"]
        rowCount = gridProperties["rowCount"]

        columnCount = gridProperties["columnCount"]

        result = (rowCount, columnCount)

        self.size_cache[key] = result
        return result

    def _get_sheet_name_from_sheet_id(self) -> str:
        key = (self.document_id, self.sheet_id)
        if key in self.sheet_name_cache:
            return self.sheet_name_cache[key]

        s = self.client.spreadsheets().get(spreadsheetId=self.document_id).execute()
        self.request_counts[self.document_id]["read"] += 1
        sheet = self._get_sheet_by_id(s["sheets"], self.sheet_id)
        result = sheet["properties"]["title"]
        self.sheet_name_cache[key] = result
        return result

    def read_sheet(self) -> List[List[str]]:
        key = f"doc:{self.document_id} sheet:{self.sheet_id}"
        cache = self.cache
        if key in cache:
            sheet = cache[key]
        else:
            assert "`" not in self.sheet_id
            sheet_name = self._get_sheet_name_from_sheet_id()
            result = (
                self.client.spreadsheets()
                .values()
                .get(spreadsheetId=self.document_id, range=f"{sheet_name}")
                .execute()
            )
            self.request_counts[self.document_id]["read"] += 1
            sheet = result["values"]
            # maybe do ragged processing here?
            cache[key] = sheet
        return sheet

    def read_row(self, row_index, start_column, end_column):
        sheet = self.read_sheet()
        row = sheet[row_index][start_column: end_column + 1]
        row.extend([""] * (end_column - len(row)))
        assert len(row) == end_column - start_column
        return row

    def read_column(self, column_index, start_row, end_row):
        sheet = self.read_sheet()
        rows = sheet[start_row:end_row]
        column = []
        for i in range(end_row - start_row):
            if i < len(rows):
                if column_index < len(rows[i]):
                    value = rows[i][column_index]
                else:
                    value = ""
            else:
                value = ""
            column.append(value)
        return column

    def write_column(
        self, column_index, start_row, values,
    ):
        assert isinstance(
            values, List
        ), f"values is of type {type(values)} but expected List"
        sheet_name = self._get_sheet_name_from_sheet_id()
        end_row = start_row + len(values)
        cell_range = "{sheet_name}!{column}{start_row}:{column}{end_row}".format(
            sheet_name=sheet_name,
            start_row=start_row + 1,
            column=colnum_string(column_index + 1),
            end_row=end_row + 1,
        )

        value_range_body = {
            "range": cell_range,  # maybe needs to drop {sheet name} prefix?
            "majorDimension": "COLUMNS",
            "values": [values],
        }

        request = (
            self.client.spreadsheets()
            .values()
            .update(
                spreadsheetId=self.document_id,
                range=cell_range,
                valueInputOption="RAW",
                body=value_range_body,
            )
        )
        self.request_counts[self.document_id]["write"] += 1

        response = request.execute()

    # Commenting out since don't have google drive api access to master file
    # def get_owners(self, self.document_id: str) -> List[Tuple[str, str]]:
    #     request = self.drive_service.files().get(fileId=self.document_id, fields="owners")
    #     response = request.execute()
    #     return [(o["displayName"], o["emailAddress"]) for o in response["owners"]]
