"""BRENDA API client
"""
# Make it organism specific? Not sure if necessary and doesn't seem to be
# well curated, at least not for mouse

from zeep import Client
import hashlib


class BRENDA_api():
    """
    API to access BREANDA entries
    """

    def __init__(
        self,
        email: str = None,
        password: str = None,
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
            input(
                "Enter your BRENDA password: "
                ).encode("utf-8")).hexdigest()

    def intialize_client(self):
        wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
        self.client = Client(wsdl)

    def get_substrates(self, ec_number: str, organism: str = '') -> list:
        """
        Get a list of substrates for a given EC number
        """
        substrates = []
        parameters = (
            self.email,
            self.password,
            f"ecNumber*{ec_number}",
            f"organism*{organism}",
            "naturalSubstrate*",
            "naturalReactionPartners*",
            "ligandStructureId*",
            )
        resultString = self.client.service.getNaturalSubstrate(*parameters)
        for entry in resultString:
            [
                substrates.append(entry)
                for entry in entry.naturalSubstrate.split(" + ")
                if entry not in ['more', '?']
                ]
        substrates = list(set(substrates))
        return substrates

    def get_products(self, ec_number: str, organism: str = '') -> list:
        """
        Get a list of products for a given EC number
        """
        products = []
        parameters = (
            self.email,
            self.password,
            f"ecNumber*{ec_number}",
            f"organism*{organism}",
            "naturalProduct*",
            "naturalReactionPartners*",
            "ligandStructureId*",
            )
        resultString = self.client.service.getNaturalProduct(*parameters)
        for entry in resultString:
            [
                products.append(entry)
                for entry in entry.naturalProduct.split(" + ")
                if entry not in ['more', '?']
                ]
        products = list(set(products))
        return products
