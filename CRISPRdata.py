#!/usr/bin/env python
# A mixture of data from sx, lcy, sj, ljh, lier

import re
import datasets
import re._compiler
from typing import Callable
import os
import torch
import torch.nn.functional as F

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

class CRISPRdataConfig(datasets.BuilderConfig):
    def __init__(self, trans_func: Callable | None = None, ref_filter: Callable | None = None, cut_filter: Callable | None = None, author_filter: Callable | None = None, file_filter: Callable | None = None, test_ratio: float = 0.05, validation_ratio: float = 0.05, seed: int = 63036, features: datasets.Features = None, **kwargs):
        """BuilderConfig for CRISPRdata.
        Args:
        trans_func: *function*, transform function applied after filter.
        ref_filter: *function*, ref_filter(ref1, ref2) -> bool.
        cut_filter: *function*, cut_filter(cut1, cut2, ref1[optional], ref2[optional]) -> bool.
        author_filter: *function*, author_filter(author, ref1[optional], ref2[optional], cut1[optional], cut2[optional]) -> bool.
        file_filter: *function*, file_filter(file, ref1[optional], ref2[optional], cut1[optional], cut2[optional], author[optional]) -> bool.
        test_ratio: *float*, the ratio of data for test.
        validation_ratio: *float*, the ratio of data for validation.
        seed: *int*, the random seed.
        **kwargs: keyword arguments forwarded to super.
        """
        super().__init__(**kwargs)
        self.trans_func = trans_func
        self.ref_filter = ref_filter
        self.cut_filter = cut_filter
        self.author_filter = author_filter
        self.file_filter = file_filter
        self.test_ratio = test_ratio
        self.validation_ratio = validation_ratio
        self.seed = seed
        self.features = features

class CRISPRdata(datasets.GeneratorBasedBuilder):
    class CRISPRinnderData:
        def __init__(self, ref1len: int = 127, ref2len: int = 127) -> None:
            self.ref1len = ref1len
            self.ref2len = ref2len
            self.features_inDelphi = datasets.Features({
                'ref': datasets.Value('string'),
                'cut': datasets.Sequence(datasets.Value('uint16')),
                'mh_gt_pos': datasets.Sequence(datasets.Value('uint16')),
                'mh_del_len': datasets.Sequence(datasets.Value('uint16')),
                'mh_mh_len': datasets.Sequence(datasets.Value('uint16')),
                'mh_gc_frac': datasets.Sequence(datasets.Value('float32')),
                'mh_count': datasets.Sequence(datasets.Value('uint64')),
                'mhless_count': datasets.Sequence(datasets.Value('uint64')),
                'insert_1bp': datasets.Sequence(datasets.Value('uint64')),
            })
            # given offset in range(-ref2len, ref1len + 1)
            # row_idx_start = max(-offset, 0)
            # row_idx_end = min(ref1len-offset, ref2len)
            # col_idx_start = max(offset, 0)
            # col_idx_end = min(offset+ref2len, ref1len)
            # diag_indices = (
            #   tensor((ref2len+1)*(ref1len+1),),
            #   tensor((ref2len+1)*(ref1len+1),)
            # )
            self.alphacode = torch.frombuffer("ACGTN".encode(), dtype=torch.uint8)
            self.base_idx = [self.alphacode % 5]
            for i in range(1, 2):
                self.base_idx.append(
                    (self.base_idx[-1] + ((self.alphacode % 5) * len(self.alphacode) ** i)[:, None]).flatten()
                )
            cutoff = 0
            for bi in self.base_idx:
                bi += cutoff
                cutoff += len(bi)
            self.base_idx = torch.cat(self.base_idx)
            self.diag_indices = (
                torch.cat([torch.arange(max(-offset, 0), min(self.ref1len - offset, self.ref2len) + 1) for offset in range(-self.ref2len, self.ref1len + 1)]),
                torch.cat([torch.arange(max(offset, 0), min(offset + self.ref2len, self.ref1len) + 1) for offset in range(-self.ref2len, self.ref1len + 1)])
            )

        # auxiliary methods
        def get_observations(self, ref1, ref2, cut, INSERT_LIMIT=None):
            assert len(ref1) == self.ref1len and len(ref2) == self.ref2len, "reference length does not fit"
            observations = torch.zeros(len(ref2) + 1, len(ref1) + 1, dtype=torch.uint64)
            if INSERT_LIMIT is not None:
                base = len(self.alphacode)
                total = (base ** (INSERT_LIMIT + 1) - base) // (base - 1)
                assert total <= len(self.base_idx), "INSERT_LIMIT is beyonded"
                insert_counts = torch.zeros(total)
                base_vec = base ** torch.arange(INSERT_LIMIT)
                base_cutoff = (base ** (torch.arange(INSERT_LIMIT) + 1) - base) // (base - 1)
                cut1, cut2 = cut['cut1'], cut['cut2']
            for author in cut['authors']:
                for file in author['files']:
                    observations[file['ref2_start'], file['ref1_end']] += torch.tensor(file['count'])
                    if INSERT_LIMIT is not None:
                        for ref1_end, ref2_start, random_insert, count in zip(file['ref1_end'], file['ref2_start'], file['random_insert'], file['count']):
                            if ref1_end < cut1 or ref2_start > cut2:
                                continue
                            insert = ref1[cut1:ref1_end] + random_insert + ref2[ref2_start:cut2]
                            if len(insert) > 0 and len(insert) <= INSERT_LIMIT:
                                insert_counts[base_cutoff[len(insert) - 1] + base_vec.inner(torch.frombuffer(insert.encode(), dtype=torch.uint8))] += count
            if INSERT_LIMIT is not None:
                return observations, insert_counts[self.base_idx[:len(insert_counts)]]
            return observations

        def num2micro_homology(self, ref1num, ref2num, cut1, cut2, observations, output_mh_merge=False):
            # OUTPUT: micro-homology matrix shape=(ref2len+1, ref1len+1) dtype=int64
            assert len(ref1num) == self.ref1len and len(ref2num) == self.ref2len, "reference length does not fit"
            indices_num = len(self.diag_indices[0])
            mh_matrix = F.pad((ref1num[:cut1].view(1, -1) == ref2num[cut2:].view(-1, 1)).to(torch.int64), pad=(0, len(ref1num) - cut1 + 1, cut2, 1), value=0)
            rep_num = torch.cat((
                torch.tensor([-1]),
                torch.where(mh_matrix[self.diag_indices].diff())[0],
                torch.tensor([indices_num - 1])
            )).diff()
            rep_val = rep_num.clone()
            rep_val[0::2] = 0
            rep_num[1::2] = rep_num[1::2] + 1
            rep_num[2::2] = rep_num[2::2] - 1
            rep_val = rep_val.repeat_interleave(rep_num)
            mh_matrix[self.diag_indices] = rep_val
            if output_mh_merge:
                gt_poss = (torch.arange(-cut2, len(ref2num) - cut2 + 1) + cut1)[:, None].expand(-1, len(ref1num) + 1)[self.diag_indices]
                del_lens = (torch.arange(cut1, cut1 - len(ref1num) - 1, -1)[None, :] + torch.arange(-cut2, len(ref2num) - cut2 + 1)[:, None])[self.diag_indices]
                mh_idxs = rep_num.cumsum(dim=0)[1::2] - 1
                mask = rep_val == 0
                del_lens = torch.cat([del_lens[mask], del_lens[mh_idxs]])
                gt_poss = torch.cat([gt_poss[mask], gt_poss[mh_idxs]])
                mh_lens = torch.cat([
                    mh_matrix[self.diag_indices][mask],
                    mh_matrix[self.diag_indices][mh_idxs]
                ])
                counts = torch.zeros(len(rep_num) // 2, dtype=torch.uint64)
                counts = counts.scatter_add(
                    dim = 0,
                    index = torch.arange(len(rep_num) // 2).repeat_interleave(rep_num[1::2]),
                    src = observations[self.diag_indices][~mask]
                )
                counts = torch.cat([observations[self.diag_indices][mask], counts])
                return del_lens, mh_lens, gt_poss, counts
            return mh_matrix

        def DNA2num(self, seq):
            # INPUT: string
            # OUTPUT: tensor dtype=uint8 shape=(len(seq),)
            seq = torch.frombuffer(seq.encode(), dtype=torch.uint8)
            cover = torch.full((len(seq),), False)
            for i in range(len(self.alphacode)):
                mask = seq == self.alphacode[i]
                seq[mask] = i
                cover = cover | mask
            assert torch.all(cover), "find unknow alphabet"
            return seq

        # trans_funcs
        def inDelphi_trans_func(self, examples, DELLEN_LIMIT=60):
            def gc_content(seq):
                return (seq.count("G") + seq.count("C")) / len(seq)
            refs, cut_list, mh_del_lenss, mh_mh_lenss, mh_gt_posss, mh_gc_fracss, mh_countss, mhless_countss, insert_1bpss = [], [], [], [], [], [], [], [], []
            for ref1, ref2, cuts in zip(examples['ref1'], examples['ref2'], examples['cuts']):
                for cut in cuts:
                    cut1, cut2 = cut['cut1'], cut['cut2']
                    refs.append(ref1[:cut1] + ref2[cut2:])
                    cut_list.append(cut1)
                    observations, insert_counts = self.get_observations(ref1, ref2, cut, INSERT_LIMIT=1)
                    insert_1bpss.append(insert_counts)
                    del_lens, mh_lens, gt_poss, counts = self.num2micro_homology(self.DNA2num(ref1), self.DNA2num(ref2), cut1, cut2, observations, output_mh_merge=True)
                    breakpoint()
                    mask = del_lens > 0 and del_lens < DELLEN_LIMIT and gt_poss >= cut1 and gt_poss - del_lens <= cut1
                    del_lens, mh_lens, gt_poss, counts = del_lens[mask], mh_lens[mask], gt_poss[mask], counts[mask]
                    mask = mh_lens > 0
                    mh_del_lenss.append(del_lens[mask])
                    mh_gt_posss.append(gt_poss[mask])
                    mh_gc_fracss.append(torch.tensor([gc_content(refs[-1][gt_pos - del_len:gt_pos]) for del_len, gt_pos in zip(mh_del_lenss[-1], mh_gt_posss[-1])]))
                    mh_countss.append(counts[mask])
                    mh_mh_lenss.append(mh_lens[mask])
                    mhless_counts = torch.zeros(DELLEN_LIMIT - 1, dtype=torch.uint64)
                    mhless_counts = mhless_counts.scatter_add(dim = 0, index=del_lens[~mask] - 1, src=counts[~mask])
                    mhless_countss.append(mhless_counts)
            return {
                'ref': refs,
                'cut': cut_list,
                'mh_gt_pos': mh_gt_posss,
                'mh_del_len': mh_del_lenss,
                'mh_mh_len': mh_mh_lenss,
                'mh_gc_frac': mh_gc_fracss,
                'mh_count': mh_countss,
                'mhless_count': mhless_countss,
                'insert_1bp': insert_1bpss
            }

    inner_data = CRISPRinnderData()

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
        CRISPRdataConfig(
            trans_func = inner_data.inDelphi_trans_func,
            author_filter=lambda author, ref1, ref2, cut1, cut2: author == "SX",
            file_filter=lambda file, ref1, ref2, cut1, cut2, author: bool(re.search("^(A2-|A7-|D2-)", file)),
            features = inner_data.features_inDelphi,
            name="SX_spcas9_inDelphi",
            version=VERSION,
            description="Data of spcas9 protein of sx and lcy for inDelphi training",
        ),
        CRISPRdataConfig(
            trans_func = inner_data.inDelphi_trans_func,
            author_filter=lambda author, ref1, ref2, cut1, cut2: author == "SX",
            file_filter=lambda file, ref1, ref2, cut1, cut2, author: bool(re.search("^(X-|x-|B2-|36t-)", file)),
            features = inner_data.features_inDelphi,
            name="SX_spymac_inDelphi",
            version=VERSION,
            description="Data of spymac protein of sx and lcy for inDelphi training",
        ),
        CRISPRdataConfig(
            trans_func = inner_data.inDelphi_trans_func,
            author_filter=lambda author, ref1, ref2, cut1, cut2: author == "SX",
            file_filter=lambda file, ref1, ref2, cut1, cut2, author: bool(re.search("^(i10t-|i83-)", file)),
            features = inner_data.features_inDelphi,
            name="SX_ispymac_inDelphi",
            version=VERSION,
            description="Data of ispymac protein of sx and lcy for inDelphi training",
        ),
    ]

    # DEFAULT_CONFIG_NAME = "SX_spcas9"  # It's not mandatory to have a default configuration. Just use one if it make sense.

    def _info(self):
        return datasets.DatasetInfo(
            # This is the description that will appear on the datasets page.
            description="""\
                This dataset is used to train a DL model predicting editing results of CRISPR.
            """,
            # This defines the different columns of the dataset and their types
            # features=
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

    def _split_generators(self, dl_manager):
        # TODO: This method is tasked with downloading/extracting the data and defining the splits depending on the configuration
        # If several configurations are possible (listed in BUILDER_CONFIGS), the configuration selected by the user is in self.config.name

        # dl_manager is a datasets.download.DownloadManager that can be used to download and extract URLS
        # It can accept any type or nested list/dict and will give back the same structure with the url replaced with path to local files.
        # By default the archives will be extracted and a path to a cached folder where they are extracted is returned instead of the archive
        hf_endpoint = os.environ.get("HF_ENDPOINT", "https://" + "huggingface.co")
        downloaded_files = dl_manager.download(f"{hf_endpoint}/datasets/ljw20180420/CRISPRdata/resolve/main/dataset.json.gz")
        ds = datasets.load_dataset('json', data_files=downloaded_files)
        ds = ds.map(self.filter_refs, batched=True)
        if self.config.trans_func is not None:
            ds = ds.map(self.config.trans_func, batched=True)
        ds = self.split_train_valid_test(ds)
        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "dataset": ds['train'],
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.VALIDATION,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "dataset": ds['validation'],
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.TEST,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "dataset": ds['test'],
                },
            )
        ]

    # method parameters are unpacked from `gen_kwargs` as given in `_split_generators`
    def _generate_examples(self, dataset):
        # TODO: This method handles input defined in _split_generators to yield (key, example) tuples from the dataset.
        # The `key` is for legacy reasons (tfds) and is not important in itself, but must be unique for each example.
        for id, example in enumerate(dataset):
            yield id, example

    def split_train_valid_test(self, ds):
        # Divide ds's train split to validation and test splits.
        ds = ds['train'].train_test_split(test_size=self.config.test_ratio + self.config.validation_ratio, shuffle=True, seed=self.config.seed) 
        ds_valid_test = ds['test'].train_test_split(test_size=self.config.test_ratio / (self.config.test_ratio + self.config.validation_ratio), shuffle=False)
        ds['validation'] = ds_valid_test.pop('train')
        ds['test'] = ds_valid_test.pop('test')
        return ds
    
    def filter_refs(self, examples):
        ref1s, ref2s, cutss = [], [], []
        for ref1, ref2, cuts in zip(examples['ref1'], examples['ref2'], examples['cuts']):
            if self.config.ref_filter is None or self.config.ref_filter(ref1, ref2):
                if self.config.cut_filter is not None or self.config.author_filter is not None or self.config.file_filter is not None:
                    cuts = self.filter_cuts(cuts, ref1, ref2)
                if cuts:
                    ref1s.append(ref1)
                    ref2s.append(ref2)
                    cutss.append(cuts)
        return {
            "ref1": ref1s,
            "ref2": ref2s,
            "cuts": cutss
        }

    def filter_cuts(self, cuts, ref1, ref2):
        new_cuts = []
        for cut in cuts:
            if self.config.cut_filter is None or self.config.cut_filter(cut["cut1"], cut["cut2"], ref1, ref2):
                if self.config.author_filter is not None or self.config.file_filter is not None:
                    cut["authors"] = self.filter_authors(cut["authors"], ref1, ref2, cut["cut1"], cut["cut2"])
                if cut["authors"]:
                    new_cuts.append(cut)
        return new_cuts

    def filter_authors(self, authors, ref1, ref2, cut1, cut2):
        new_authors = []
        for author in authors:
            if self.config.author_filter is None or self.config.author_filter(author["author"], ref1, ref2, cut1, cut2):
                if self.config.file_filter is not None:
                    author["files"] = self.filter_files(author["files"], ref1, ref2, cut1, cut2, author["author"])
                if author["files"]:
                    new_authors.append(author)
        return new_authors
    
    def filter_files(self, files, ref1, ref2, cut1, cut2, author):
        return [file for file in files if self.config.file_filter(file["file"], ref1, ref2, cut1, cut2, author)]


if __name__ == "__main__":
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
    
    # from datasets.commands.test import TestCommand
    # test_command = TestCommand(
    #     dataset="./CRISPRdata.py",
    #     name=None,
    #     cache_dir=None,
    #     data_dir=None,
    #     all_configs=True,
    #     save_infos=True,
    #     ignore_verifications=False,
    #     force_redownload=False,
    #     clear_cache=False,
    #     num_proc=None,
    #     trust_remote_code=True
    # )
    # test_command.run()

    # api.upload_file(
    #     path_or_fileobj="./README.md",
    #     path_in_repo="README.md",
    #     repo_id="ljw20180420/CRISPRdata",
    #     repo_type="dataset",
    # )
