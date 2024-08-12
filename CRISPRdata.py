#!/usr/bin/env python
# A mixture of data from sx, lcy, sj, ljh, lier

import re
import datasets

# TODO: Add BibTeX citation
# Find for instance the citation on arxiv or on the dataset repo/website
# _CITATION = """\
# @InProceedings{huggingface:dataset,
# title = {A great new dataset},
# author={huggingface, Inc.
# },
# year={2020}
# }
# """

# TODO: Add a link to an official homepage for the dataset here
# _HOMEPAGE = ""

# TODO: Add the licence for the dataset here if you can find it
# _LICENSE = ""

# TODO: Add link to the official dataset URLs here
# The HuggingFace Datasets library doesn't host the datasets but only points to the original files.
# This can be an arbitrary nested dict/list of URLs (see below in `_split_generators` method)
# _URLS = {
#     "first_domain": "https://huggingface.co/great-new-dataset-first_domain.zip",
#     "second_domain": "https://huggingface.co/great-new-dataset-second_domain.zip",
# }

def spcas9_sxlcy_filter(examples):
    return [
        author == "SX" and re.search("^(A2-|A7-|D2-)", file_name) for author, file_name in zip(examples['author'], examples['file_name'])
    ]

def spymac_sxlcy_filter(examples):
    return [
        author == "SX" and re.search("^(X-|x-|B2-|36t-)", file_name) for author, file_name in zip(examples['author'], examples['file_name'])
    ]

def ispymac_sxlcy_filter(examples):
    return [
        author == "SX" and re.search("^(i10t-|i83-)", file_name) for author, file_name in zip(examples['author'], examples['file_name'])
    ]

class CRISPRdataConfig(datasets.BuilderConfig):
    def __init__(self, filter_func, test_valid_ratio=0.1, seed=63036, **kwargs):
        """BuilderConfig for CRISPRdata.

        Args:
        filter_func: *function*, filter function for configurations.
        **kwargs: keyword arguments forwarded to super.
        """
        super().__init__(**kwargs)
        self.filter_func = filter_func
        self.test_valid_ratio=test_valid_ratio
        self.seed=seed

class CRISPRdata(datasets.ArrowBasedBuilder):
    # This is an example of a dataset with multiple configurations.
    # If you don't want/need to define several sub-sets in your dataset,
    # just remove the BUILDER_CONFIG_CLASS and the BUILDER_CONFIGS attributes.

    # If you need to make complex sub-parts in the datasets with configurable options
    # You can create your own builder configuration class to store attribute, inheriting from datasets.BuilderConfig
    # BUILDER_CONFIG_CLASS = MyBuilderConfig

    # You will be able to load one or the other configurations in the following list with
    # data = datasets.load_dataset('path_to_CRISPRdata', 'config_name')

    VERSION = datasets.Version("1.0.0")

    BUILDER_CONFIG_CLASS = CRISPRdataConfig

    BUILDER_CONFIGS = [
        CRISPRdataConfig(spcas9_sxlcy_filter, name="spcas9_sxlcy", version=VERSION, description="Data of spcas9 protein of sx and lcy"),
        CRISPRdataConfig(spymac_sxlcy_filter, name="spymac_sxlcy", version=VERSION, description="Data of spymac protein of sx and lcy"),
        CRISPRdataConfig(ispymac_sxlcy_filter, name="ispymac_sxlcy", version=VERSION, description="Data of ispymac protein of sx and lcy"),
    ]

    # DEFAULT_CONFIG_NAME = "spcas9_sxlcy"  # It's not mandatory to have a default configuration. Just use one if it make sense.

    def _info(self):
        return datasets.DatasetInfo(
            # This is the description that will appear on the datasets page.
            description="""\
                This dataset is used to train a DL model predicting editing results of CRISPR.
            """,
            # This defines the different columns of the dataset and their types
            features=datasets.Features({
                'ref1': datasets.Value('string'),
                'ref2': datasets.Value('string'),
                'cut1': datasets.Sequence(datasets.Value('uint16')),
                'cut2': datasets.Sequence(datasets.Value('uint16')),
                'ref1_end': datasets.Sequence(datasets.Sequence(datasets.Value('uint16'))),
                'ref2_start': datasets.Sequence(datasets.Sequence(datasets.Value('uint16'))),
                'random_insert': datasets.Sequence(datasets.Sequence(datasets.Value('string'))),
                'count': datasets.Sequence(datasets.Sequence(datasets.Value('uint64'))),
                'author': datasets.Value('string'),
                'file_name': datasets.Value('string')
            })
            # If there's a common (input, target) tuple from the features, uncomment supervised_keys line below and
            # specify them. They'll be used if as_supervised=True in builder.as_dataset.
            # supervised_keys=("sentence", "label"),
            # Homepage of the dataset for documentation
            # homepage=_HOMEPAGE,
            # License for the dataset if available
            # license=_LICENSE,
            # Citation for the dataset
            # citation=_CITATION,
        )

    def split_train_valid_test(self, ds):
        # Divide ds's train split to valid and test splits. Both has proportion test_valid_ratio.
        ds = ds['train'].train_test_split(2 * self.config.test_valid_ratio, shuffle=True, seed=self.config.seed) 
        ds_valid_test = ds['test'].train_test_split(test_size=0.5, shuffle=False)
        ds['validation'] = ds_valid_test.pop('train')
        ds['test'] = ds_valid_test.pop('test')
        return ds

    def _split_generators(self, dl_manager):
        # TODO: This method is tasked with downloading/extracting the data and defining the splits depending on the configuration
        # If several configurations are possible (listed in BUILDER_CONFIGS), the configuration selected by the user is in self.config.name

        # dl_manager is a datasets.download.DownloadManager that can be used to download and extract URLS
        # It can accept any type or nested list/dict and will give back the same structure with the url replaced with path to local files.
        # By default the archives will be extracted and a path to a cached folder where they are extracted is returned instead of the archive
        ds = datasets.load_dataset('json', data_files="dataset.json.gz")
        ds = ds.filter(self.config.filter_func, batched=True)
        ds = self.split_train_valid_test(ds)
        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "dataset": ds['train']
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.VALIDATION,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "dataset": ds['validation']
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.TEST,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "dataset": ds['test']
                },
            ),
        ]

    # method parameters are unpacked from `gen_kwargs` as given in `_split_generators`
    def _generate_tables(self, dataset):
        # TODO: This method handles input defined in _split_generators to yield (key, example) tuples from the dataset.
        # The `key` is for legacy reasons (tfds) and is not important in itself, but must be unique for each example.
        for batch in dataset.iter(batch_size=100):
            yield 0, batch

if __name__ == "__main__":
    from datasets.commands.test import TestCommand
    test_command = TestCommand(
        dataset="./CRISPRdata.py",
        name=None,
        cache_dir=None,
        data_dir=None,
        all_configs=True,
        save_infos=True,
        ignore_verifications=False,
        force_redownload=False,
        clear_cache=False,
        num_proc=None,
        trust_remote_code=True
    )
    test_command.run()

    from huggingface_hub import HfApi
    api = HfApi()
    api.upload_file(
        path_or_fileobj="./CRISPRdata.py",
        path_in_repo="CRISPRdata.py",
        repo_id="ljw20180420/CRISPRdata",
        repo_type="dataset",
    )
    api.upload_file(
        path_or_fileobj="./dataset.json.gz",
        path_in_repo="dataset.json.gz",
        repo_id="ljw20180420/CRISPRdata",
        repo_type="dataset",
    )
    api.upload_file(
        path_or_fileobj="./README.md",
        path_in_repo="README.md",
        repo_id="ljw20180420/CRISPRdata",
        repo_type="dataset",
    )
