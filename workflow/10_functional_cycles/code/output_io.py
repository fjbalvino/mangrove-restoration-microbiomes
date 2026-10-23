"""Write complete binary outputs atomically and verify their stored contents."""
from pathlib import Path
from io import BytesIO
import hashlib
import os
import tempfile
import numpy as np


def write_binary(path, payload):
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.'+path.name+'.', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        if hashlib.sha256(path.read_bytes()).digest() != hashlib.sha256(payload).digest():
            raise OSError('Written output does not match its complete in-memory payload: '+str(path))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def save_npz(path, **arrays):
    stream = BytesIO()
    np.savez_compressed(stream, **arrays)
    write_binary(path, stream.getvalue())
