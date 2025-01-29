from distutils.command.build_ext import build_ext as build_ext_orig
from setuptools import Extension
from setuptools import setup
#  https://stackoverflow.com/questions/4529555/building-a-ctypes-based-c-library-with-distutils


class CTypesExtension(Extension):
    pass


class build_ext(build_ext_orig):
    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypesExtension)
        return super().build_extension(ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return super().get_export_symbols(ext)

    def get_ext_filename(self, ext_name):
        if self._ctypes:
            return ext_name + ".so"
        return super().get_ext_filename(ext_name)


setup(
    ext_modules=[
        CTypesExtension(
            "ngsTroubleFinder._pyPaCBAM",
            ["src/ngsTroubleFinder/pyPaCBAM.c"],
            libraries=["hts"],
        ),
    ],
    cmdclass={"build_ext": build_ext},
)
