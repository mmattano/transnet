"""SQLite3 database functions
"""

__all__ = [
    "create_connection",
    "execute_query",
    "create_table_query",
    "list_tables_in_sqlite3_db",
]

import sqlite3
from sqlite3 import Error
from pypika import Query


def create_connection(path):
    """
    Create a database connection to the SQLite database specified by db_file.
    If no file exists, it will be created.

    Parameters
    ----------
    path : str
        Either path to an existing database or path to a new database to be
        created.

    Returns
    -------
    connection : sqlite3.Connection

    Raises
    ------
    Error : sqlite3.Error
    """

    connection = None
    try:
        connection = sqlite3.connect(path)
        print("Connection to SQLite DB successful")
    except Error as e:
        print(f"The error '{e}' occurred")

    return connection


def execute_query(connection, query, output=False):
    """
    Execute a query on the database.

    Parameters
    ----------
    connection : sqlite3.Connection
        A connection to the database, initialized e.g. with create_connection()
    query : str
        The query to execute written in SQL format.

    Raises
    ------
    Error : sqlite3.Error
    """

    cursor = connection.cursor()
    try:
        cursor.execute(query)
        connection.commit()
        print("Query executed successfully")
        if output:
            return cursor.fetchall()
    except Error as e:
        print(f"The error '{e}' occurred")


def create_table_query(table, columns_list):
    """
    Create a new table in the database.

    Parameters
    ----------
    table : pypika.Table
        A table element with the name of the table to be created.
    columns_list : list
        List of pypika.Column elements.
    """

    return Query.create_table(table).columns(*columns_list).get_sql()


def list_tables_in_sqlite3_db(connection):
    """
    List all tables in the database.

    Parameters
    ----------
    connection : sqlite3.Connection
        A connection to the database, initialized e.g. with
        create_connection().

    Returns
    -------
    tables : list
        List of tables in the database.
    """

    cursor = connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table';"
    )
    tables = [v[0] for v in cursor.fetchall() if v[0] != "sqlite_sequence"]
    cursor.close()
    return tables
