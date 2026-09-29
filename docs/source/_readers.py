"""Read the jupytext walkthroughs, with a kernel that exists.

``notebooks/*.py`` carry no kernel spec, so myst-nb falls back to whichever
kernel is registered as ``python3`` -- which on a machine with several conda
environments is usually the wrong one, and sometimes one that has been
deleted. This names a kernel of our own and installs it if it is missing, so
``make docs`` works from a fresh checkout.
"""

KERNEL = "transnet"


def _ensure_kernel() -> str:
    from jupyter_client.kernelspec import KernelSpecManager

    if KERNEL not in KernelSpecManager().find_kernel_specs():
        from ipykernel.kernelspec import install

        install(user=True, kernel_name=KERNEL, display_name="TransNet")
    return KERNEL


def read(text: str):
    """Parse a jupytext light-format script into a notebook."""
    import jupytext

    notebook = jupytext.reads(text, fmt="py:light")
    notebook.metadata["kernelspec"] = {
        "name": _ensure_kernel(), "display_name": "TransNet", "language": "python",
    }
    return notebook
