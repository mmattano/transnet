"""BRENDA API client
"""
# Make it organism specific? Not sure if necessary and doesn't seem to be
# well curated, at least not for mouse

from zeep import Client
from zeep.helpers import serialize_object
from zeep.exceptions import TransportError
import hashlib
import time


class BRENDA_api:
    """
    API to access BREANDA entries
    """

    def __init__(
        self, email: str = None, password: str = None,
    ):
        self.email = email
        self.password = password
        if self.email is None and self.password is None:
            self.sign_in()
        elif self.email is None:
            self.change_email()
        elif self.password is None:
            self.change_password()
        self.intialize_client()

    def sign_in(self):
        self.change_email()
        self.change_password()

    def change_email(self):
        self.email = input("Enter your BRENDA account's email: ")

    def change_password(self):
        self.password = hashlib.sha256(
            input("Enter your BRENDA password: ").encode("utf-8")
        ).hexdigest()

    def intialize_client(self):
        wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
        self.client = Client(wsdl)