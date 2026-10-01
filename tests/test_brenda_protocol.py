"""BRENDA's SOAP calling convention.

BRENDA does not take ordinary named parameters. Every field after the
credentials must be a positional ``"field*value"`` string, and every field the
method declares must be present. Getting this wrong is not an error -- BRENDA
returns an empty result, indistinguishable from "this enzyme has no data".
That is how 3,334 queries came back empty without a single exception.
"""

import pytest

pytest.importorskip("zeep", reason="pip install zeep")

from transnet.api.brenda import BrendaClient          # noqa: E402


class FakeOperation:
    def __init__(self, names):
        self.input = type("Body", (), {
            "body": type("Elem", (), {
                "type": type("CT", (), {
                    "elements": [(n, None) for n in names]
                })()
            })()
        })()


def _client_with(monkeypatch, names, capture):
    """A BrendaClient whose SOAP service records the arguments it is given."""
    client = BrendaClient(email="user@example.com", password="secret")

    class Service:
        def getInhibitors(self, *args):
            capture["args"] = args
            return []

    class Port:
        binding = type("B", (), {"_operations": {"getInhibitors": FakeOperation(names)}})()

    fake = type("Client", (), {
        "service": Service(),
        "wsdl": type("W", (), {
            "services": {"s": type("S", (), {"ports": {"p": Port()}})()}
        })(),
    })()
    monkeypatch.setattr(client, "_get_client", lambda: fake)
    return client


NAMES = ["email", "password", "ecNumber", "organism",
         "inhibitor", "commentary", "ligandStructureId", "literature"]


class TestCallConvention:
    def test_fields_are_passed_as_field_star_value(self, monkeypatch):
        capture = {}
        client = _client_with(monkeypatch, NAMES, capture)
        client.get_inhibitors("1.1.1.1", "Mus musculus")

        args = capture["args"]
        assert args[2] == "ecNumber*1.1.1.1", (
            "a plain value makes BRENDA return nothing, silently"
        )
        assert args[3] == "organism*Mus musculus"

    def test_credentials_are_passed_raw(self, monkeypatch):
        capture = {}
        client = _client_with(monkeypatch, NAMES, capture)
        client.get_inhibitors("1.1.1.1", None)

        assert capture["args"][0] == "user@example.com"
        assert "*" not in capture["args"][0]
        # The password is hashed, never sent in the clear.
        assert capture["args"][1] != "secret"
        assert len(capture["args"][1]) == 64

    def test_every_declared_field_is_present(self, monkeypatch):
        capture = {}
        client = _client_with(monkeypatch, NAMES, capture)
        client.get_inhibitors("1.1.1.1", None)

        assert len(capture["args"]) == len(NAMES), (
            "BRENDA requires every declared field, empty ones as 'field*'"
        )
        for argument in capture["args"][2:]:
            assert "*" in argument

    def test_unused_filters_are_sent_empty(self, monkeypatch):
        capture = {}
        client = _client_with(monkeypatch, NAMES, capture)
        client.get_inhibitors("1.1.1.1", None)

        assert "inhibitor*" in capture["args"]
        assert "commentary*" in capture["args"]
        assert "organism*" in capture["args"], "an absent organism means no filter"

    def test_arguments_are_positional_not_keyword(self, monkeypatch):
        """BRENDA's service rejects keyword arguments."""
        captured = {}
        client = BrendaClient(email="a@b.c", password="p")

        class Service:
            def getInhibitors(self, *args, **kwargs):
                captured["args"] = args
                captured["kwargs"] = kwargs
                return []

        class Port:
            binding = type("B", (), {
                "_operations": {"getInhibitors": FakeOperation(NAMES)}
            })()

        fake = type("Client", (), {
            "service": Service(),
            "wsdl": type("W", (), {
                "services": {"s": type("S", (), {"ports": {"p": Port()}})()}
            })(),
        })()
        monkeypatch.setattr(client, "_get_client", lambda: fake)
        client.get_inhibitors("1.1.1.1", None)

        assert captured["args"], "arguments must be positional"
        assert not captured["kwargs"], "BRENDA does not accept keyword arguments"


class TestParameterDiscovery:
    def test_names_come_from_the_wsdl_per_method(self, monkeypatch):
        capture = {}
        # getSubstrate declares different fields than getInhibitors.
        substrate_names = ["email", "password", "ecNumber", "organism",
                           "substrate", "reactionPartners", "ligandStructureId"]
        client = _client_with(monkeypatch, substrate_names, capture)
        names = client._parameter_names("getInhibitors")
        assert names == substrate_names

    def test_names_are_cached(self, monkeypatch):
        capture = {}
        client = _client_with(monkeypatch, NAMES, capture)
        first = client._parameter_names("getInhibitors")
        monkeypatch.setattr(client, "_get_client",
                            lambda: (_ for _ in ()).throw(AssertionError("refetched")))
        assert client._parameter_names("getInhibitors") == first


@pytest.mark.network
def test_the_live_wsdl_still_declares_the_expected_fields():
    """If BRENDA changes its schema, this fails rather than returning empties."""
    client = BrendaClient(email="probe@example.com", password="probe")
    names = client._parameter_names("getInhibitors")
    assert names[:4] == ["email", "password", "ecNumber", "organism"]
    assert "inhibitor" in names


class TestIllegalXmlIsRecovered:
    """BRENDA commentary sometimes carries XML-illegal control bytes.

    zeep then rejected the whole response, so every activator or inhibitor
    for that EC was lost -- and retrying returns the same bytes. 5 ECs for
    mouse and 11 for human were dropped this way.
    """

    def test_control_bytes_are_stripped_before_parsing(self, monkeypatch):
        from transnet.api import brenda
        from zeep.transports import Transport

        class Response:
            def __init__(self, content):
                self.content = content
                self._content = content

        dirty = b'<?xml version="1.0"?><a>citrate\x05 inhibits</a>'
        monkeypatch.setattr(Transport, "post", lambda self, a, m, h: Response(dirty))

        response = brenda._sanitising_transport().post("http://x", b"", {})
        import lxml.etree as etree
        assert etree.fromstring(response._content).text == "citrate inhibits"

    def test_clean_responses_are_untouched(self, monkeypatch):
        from transnet.api import brenda
        from zeep.transports import Transport

        class Response:
            def __init__(self, content):
                self.content = content
                self._content = content

        clean = b'<?xml version="1.0"?><a>tab\\tand newline\\n are legal</a>'
        monkeypatch.setattr(Transport, "post", lambda self, a, m, h: Response(clean))
        response = brenda._sanitising_transport().post("http://x", b"", {})
        assert response._content == clean
