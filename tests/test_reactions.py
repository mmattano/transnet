"""Reactions layer and stoichiometry parsing.

Stoichiometry parsing has to handle both numeric coefficients and the
non-numeric ones KEGG uses ('n' for a variable-length polymer).

The tests that call ``populate(from_api=True)`` each download about 12,000
KEGG reactions, so they are marked ``network`` and deselected by default. The
same parsing is covered offline against the checked-in reaction table.
"""

import pytest
import pandas as pd
from transnet.biology.layers import Reactions
from transnet.biology.elements import Reaction


class TestReactionsParsing:
    """Test suite for reactions parsing and stoichiometry handling."""
    
    def test_reaction_creation(self):
        """Test basic reaction object creation."""
        reaction = Reaction(
            id="R00001",
            name="Test Reaction",
            equation="A + B <=> C",
            definition="Test definition",
            enzyme=["1.1.1.1"],
            substrates=["C00001", "C00002"],
            products=["C00003"],
            stoichiometry_substrates=[1.0, 2.0],
            stoichiometry_products=[1.0]
        )
        
        assert reaction.id == "R00001"
        assert reaction.name == "Test Reaction"
        assert len(reaction.substrates) == 2
        assert len(reaction.products) == 1
        assert reaction.stoichiometry_substrates == [1.0, 2.0]
        assert reaction.stoichiometry_products == [1.0]
    
    def test_reaction_with_none_stoichiometry(self):
        """Test reaction with None (variable) stoichiometry."""
        reaction = Reaction(
            id="R00002",
            name="Variable Stoichiometry Reaction",
            equation="n A <=> B",
            definition="Variable coefficient test",
            enzyme=["1.1.1.2"],
            substrates=["C00001"],
            products=["C00002"],
            stoichiometry_substrates=[None],  # Variable stoichiometry
            stoichiometry_products=[1.0]
        )
        
        assert reaction.stoichiometry_substrates == [None]
        assert reaction.stoichiometry_products == [1.0]
    
    @pytest.mark.network
    @pytest.mark.network
    def test_reactions_layer_populate(self):
        """Test populating reactions layer from API."""
        reactions_layer = Reactions()
        reactions_layer.populate(from_api=True)
        
        assert len(reactions_layer.reactions) > 0
        
        # Check that reactions have proper structure
        sample_reaction = reactions_layer.reactions[0]
        assert hasattr(sample_reaction, 'id')
        assert hasattr(sample_reaction, 'equation')
        assert hasattr(sample_reaction, 'substrates')
        assert hasattr(sample_reaction, 'products')
        assert hasattr(sample_reaction, 'stoichiometry_substrates')
        assert hasattr(sample_reaction, 'stoichiometry_products')
    
    @pytest.mark.network
    @pytest.mark.network
    def test_stoichiometry_types(self):
        """Test that stoichiometry values are proper types."""
        reactions_layer = Reactions()
        reactions_layer.populate(from_api=True)
        
        for reaction in reactions_layer.reactions[:100]:  # Test first 100
            if reaction.stoichiometry_substrates:
                for stoic in reaction.stoichiometry_substrates:
                    # Should be either float or None
                    assert stoic is None or isinstance(stoic, (int, float))
            
            if reaction.stoichiometry_products:
                for stoic in reaction.stoichiometry_products:
                    # Should be either float or None
                    assert stoic is None or isinstance(stoic, (int, float))
    
    def test_reactions_from_dataframe(self):
        """Test creating reactions layer from DataFrame."""
        # Create test DataFrame
        test_data = {
            'reaction': ['R00001', 'R00002'],  # Use 'reaction' not 'id'
            'name': ['Reaction 1', 'Reaction 2'],
            'equation': ['A + B <=> C', 'n X <=> Y'],
            'definition': ['Def 1', 'Def 2'],
            'enzyme': [['1.1.1.1'], ['2.2.2.2']],
            'substrates': [['C00001', 'C00002'], ['C00003']],
            'products': [['C00003'], ['C00004']],
            'stoichiometry_substrates': [[1.0, 1.0], [None]],
            'stoichiometry_products': [[1.0], [1.0]]
        }
        
        df = pd.DataFrame(test_data)
        
        reactions_layer = Reactions()
        reactions_layer.populate(from_api=False, df=df)
        
        assert len(reactions_layer.reactions) == 2
        assert reactions_layer.reactions[0].id == 'R00001'
        assert reactions_layer.reactions[1].stoichiometry_substrates == [None]
    
    @pytest.mark.network
    @pytest.mark.network
    def test_reaction_equation_parsing(self):
        """Test that KEGG equation parsing handles various formats."""
        reactions_layer = Reactions()
        reactions_layer.populate(from_api=True)
        
        # Check that we have reactions with different stoichiometry patterns
        has_simple = False
        has_multiple = False
        has_variable = False
        
        for reaction in reactions_layer.reactions[:1000]:
            if reaction.stoichiometry_substrates:
                # Simple 1:1 stoichiometry
                if all(s == 1.0 for s in reaction.stoichiometry_substrates if s is not None):
                    has_simple = True
                # Multiple coefficients
                if any(s is not None and s > 1 for s in reaction.stoichiometry_substrates):
                    has_multiple = True
                # Variable stoichiometry
                if any(s is None for s in reaction.stoichiometry_substrates):
                    has_variable = True
        
        # We should have at least simple reactions
        assert has_simple, "No simple stoichiometry reactions found"


class TestReactionsSerialization:
    """Test serialization and deserialization of reactions."""
    
    def test_csv_save_load_cycle(self, tmp_path):
        """Test saving and loading reactions to/from CSV."""
        # Create test reactions
        reactions_layer = Reactions()
        test_reactions = [
            Reaction(
                id="R00001",
                name="Test 1",
                equation="A + B <=> C",
                definition="Def 1",
                enzyme=["1.1.1.1"],
                substrates=["C00001", "C00002"],
                products=["C00003"],
                stoichiometry_substrates=[1.0, 2.0],
                stoichiometry_products=[1.0]
            ),
            Reaction(
                id="R00002",
                name="Test 2",
                equation="n X <=> Y",
                definition="Def 2",
                enzyme=["2.2.2.2"],
                substrates=["C00004"],
                products=["C00005"],
                stoichiometry_substrates=[None],  # Variable
                stoichiometry_products=[1.0]
            )
        ]
        reactions_layer.reactions = test_reactions
        
        # Save to CSV (mimicking build_networks.py logic)
        reactions_data = []
        for reaction in reactions_layer.reactions:
            stoich_subs = []
            if reaction.stoichiometry_substrates:
                stoich_subs = [str(x) if x is not None else "None" 
                              for x in reaction.stoichiometry_substrates]
            
            stoich_prods = []
            if reaction.stoichiometry_products:
                stoich_prods = [str(x) if x is not None else "None" 
                               for x in reaction.stoichiometry_products]
            
            reactions_data.append({
                'reaction': reaction.id,  # Use 'reaction' not 'id'
                'name': reaction.name,
                'equation': reaction.equation,
                'definition': reaction.definition,
                'enzyme': reaction.enzyme,
                'substrates': ";".join(reaction.substrates) if reaction.substrates else "",
                'products': ";".join(reaction.products) if reaction.products else "",
                'stoichiometry_substrates': ";".join(stoich_subs) if stoich_subs else "",
                'stoichiometry_products': ";".join(stoich_prods) if stoich_prods else ""
            })
        
        df = pd.DataFrame(reactions_data)
        csv_path = tmp_path / "test_reactions.csv"
        df.to_csv(csv_path, index=False)
        
        # Load back (mimicking load_reactions_from_file logic).  Converters
        # bypass pandas' NA inference, which would otherwise read the literal
        # "None" marking variable stoichiometry as NaN and lose it.
        loaded_df = pd.read_csv(csv_path, converters={
            'stoichiometry_substrates': lambda x: x,
            'stoichiometry_products': lambda x: x,
        })
        
        def parse_stoichiometry(x):
            # str() first: pandas gives this column a numeric dtype whenever
            # every row happens to hold a single bare number.
            if pd.notna(x) and str(x):
                values = []
                for y in str(x).split(';'):
                    if y in ("None", "nan"):
                        values.append(None)
                    else:
                        try:
                            values.append(float(y))
                        except ValueError:
                            values.append(None)
                return values
            return []
        
        loaded_df['stoichiometry_substrates'] = loaded_df['stoichiometry_substrates'].apply(parse_stoichiometry)
        loaded_df['stoichiometry_products'] = loaded_df['stoichiometry_products'].apply(parse_stoichiometry)
        
        # Verify loaded data
        assert len(loaded_df) == 2
        assert loaded_df.iloc[0]['reaction'] == 'R00001'  # Use 'reaction' column
        assert loaded_df.iloc[0]['stoichiometry_substrates'] == [1.0, 2.0]
        assert loaded_df.iloc[1]['stoichiometry_substrates'] == [None]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


class TestParsingOffline:
    """The same coverage, against data/master_reactions.csv -- no network."""

    @staticmethod
    def _layer():
        from pathlib import Path

        from transnet.datasets import DATA_DIR, read_master_reactions

        if not (DATA_DIR / "master_reactions.csv").exists():
            pytest.skip("data/master_reactions.csv is not present")
        frame = read_master_reactions()
        layer = Reactions()
        layer.populate(from_api=False, df=frame)
        return layer

    def test_the_table_loads_as_reaction_objects(self):
        layer = self._layer()
        assert len(layer.reactions) > 1000
        assert all(isinstance(r, Reaction) for r in layer.reactions[:50])

    def test_stoichiometry_is_float_or_none(self):
        for reaction in self._layer().reactions[:500]:
            for value in (reaction.stoichiometry_substrates or []) + \
                    (reaction.stoichiometry_products or []):
                assert value is None or isinstance(value, float), value

    def test_variable_stoichiometry_survives_the_round_trip(self):
        """KEGG writes polymer reactions as "n C00001"; those coefficients
        come back as None rather than crashing the parse."""
        layer = self._layer()
        variable = [r for r in layer.reactions
                    if None in (r.stoichiometry_substrates or [])
                    or None in (r.stoichiometry_products or [])]
        assert variable, "no variable-stoichiometry reaction found to check"
        assert all(r.substrates or r.products for r in variable[:20])
