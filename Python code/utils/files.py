import pickle
import logging
from typing import Dict, Optional

def save_session(filename: str, data: Dict) -> bool:
    try:
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
        return True
    except Exception as e:
        logging.error(f"Save session failed: {e}")
        return False

def load_session(filename: str) -> Optional[Dict]:
    try:
        with open(filename, 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        logging.error(f"Load session failed: {e}")
        return None
