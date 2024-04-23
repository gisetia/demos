# https://pypi.org/project/goenrich/

from . import obo
from . import enrich
# import goenrich.read
# import goenrich.enrich
# import goenrich.export
# import goenrich.tools

from importlib import reload
reload(obo)
reload(enrich)
