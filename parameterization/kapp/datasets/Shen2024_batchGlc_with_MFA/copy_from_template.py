# copy files from (template) to current directory
# if you're making general edits to the procedure, make them in the (template) directory
# preserves phenotype.txt and kapp_options.py since they're model- and condition-specific
import shutil, os

shutil.copytree('../(template)', '.', dirs_exist_ok=True, ignore=shutil.ignore_patterns('phenotype.txt', 'kapp_options.py'))