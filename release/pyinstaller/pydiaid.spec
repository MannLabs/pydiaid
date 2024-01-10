# -*- mode: python ; coding: utf-8 -*-

import pkgutil
import os
import sys
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT, BUNDLE, TOC
import PyInstaller.utils.hooks
import pkg_resources
import importlib.metadata
import pydiaid


##################### User definitions
exe_name = 'pydiaid_gui'
script_name = 'pydiaid_pyinstaller.py'
if sys.platform[:6] == "darwin":
	icon = '../logos/alpha_logo.icns'
else:
	icon = '../logos/alpha_logo.ico'
block_cipher = None
location = os.getcwd()
project = "pydiaid"
remove_tests = True
bundle_name = "pydiaid"
#####################


requirements = {
	req.split()[0] for req in importlib.metadata.requires(project)
}
requirements.add(project)
requirements.add("distributed")
hidden_imports = set()
datas = []
binaries = []
checked = set()
while requirements:
	requirement = requirements.pop()
	checked.add(requirement)
	if requirement in ["pywin32"]:
		continue
	try:
		module_version = importlib.metadata.version(requirement)
	except (
		importlib.metadata.PackageNotFoundError,
		ModuleNotFoundError,
		ImportError
	):
		continue
	try:
		datas_, binaries_, hidden_imports_ = PyInstaller.utils.hooks.collect_all(
			requirement,
			include_py_files=True
		)
	except ImportError:
		continue
	datas += datas_
	# binaries += binaries_
	hidden_imports_ = set(hidden_imports_)
	if "" in hidden_imports_:
		hidden_imports_.remove("")
	if None in hidden_imports_:
		hidden_imports_.remove(None)
	requirements |= hidden_imports_ - checked
	hidden_imports |= hidden_imports_

if remove_tests:
	hidden_imports = sorted(
		[h for h in hidden_imports if "tests" not in h.split(".")]
	)
else:
	hidden_imports = sorted(hidden_imports)


hidden_imports = [h for h in hidden_imports if "__pycache__" not in h]
datas = [d for d in datas if ("__pycache__" not in d[0]) and (d[1] not in [".", "Resources", "scripts"])]

a = Analysis(
	[script_name],
	pathex=[location],
	binaries=binaries,
	datas=datas,
	hiddenimports=hidden_imports,
	hookspath=[],
	runtime_hooks=[],
	excludes=[h for h in hidden_imports if "datashader" in h],
	win_no_prefer_redirects=False,
	win_private_assemblies=False,
	cipher=block_cipher,
	noarchive=False
)
pyz = PYZ(
	a.pure,
	a.zipped_data,
	cipher=block_cipher
)

if sys.platform[:5] == "linux":
	exe = EXE(
		pyz,
		a.scripts,
		a.binaries,
		a.zipfiles,
		a.datas,
		name=bundle_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		upx_exclude=[],
		icon=icon
	)
else:
	exe = EXE(
		pyz,
		a.scripts,
		# a.binaries,
		a.zipfiles,
		# a.datas,
		exclude_binaries=True,
		name=exe_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		icon=icon
	)
	coll = COLLECT(
		exe,
		a.binaries,
		# a.zipfiles,
		a.datas,
		strip=False,
		upx=True,
		upx_exclude=[],
		name=exe_name
	)
	import shutil
	import sklearn.neighbors._partition_nodes
	import sklearn.utils._typedefs
	if sys.platform[:6] == "darwin":
		import cmath
		shutil.copyfile(
			cmath.__file__,
			f"dist/{exe_name}/{os.path.basename(cmath.__file__)}"
		)
	new_location = os.path.join(
		"dist",
		exe_name,
		"sklearn",
		"neighbors",
		os.path.basename(sklearn.neighbors._partition_nodes.__file__),
	)
	shutil.copyfile(
		sklearn.neighbors._partition_nodes.__file__,
		new_location
	)
	new_location = os.path.join(
		"dist",
		exe_name,
		"sklearn",
		"utils",
		os.path.basename(sklearn.utils._typedefs.__file__),
	)
	shutil.copyfile(
		sklearn.utils._typedefs.__file__,
		new_location
	)
