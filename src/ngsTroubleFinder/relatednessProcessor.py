class RelatednessProcessor:
    """
    Class that implements all the functions related to relatedness
    """

    def __init__(self):
        self.snpSet = {}

    def _createSNPSets(self, pileup):
        """
        Divides the autosomals variants of a sample into three different classes(ref, het, alt)
        Args:
            pileup: pd.Dataframe
                pileup in PyPaCBAM format

        Returns: Dictionary{str:set(str,int)}"
            A dictionary containing the set in form of (chromosome, position)
            of reference, heterozygous and alternative snps

        """
        pileup = pileup[pileup["cov"] >= 20]
        pileup = pileup[pileup["chr"] != "chrX"]
        pileup = pileup[pileup["chr"] != "chrY"]

        ref = pileup[pileup["genotype"] == "0/0"]
        ref = set(zip(ref["chr"], ref["pos"]))

        het = pileup[pileup["genotype"] == "0/1"]
        het = set(zip(het["chr"], het["pos"]))

        alt = pileup[pileup["genotype"] == "1/1"]
        alt = set(zip(alt["chr"], alt["pos"]))

        return {"ref": ref, "het": het, "alt": alt, "total": ref | het | alt}

    def _computeRelatedness(self, sample1, sample2):
        """
        Compute the relatedness between two samples
        Args:
            sample1: samples with a pileup
            sample2: samples with a pileup

        Returns: Tuple(float, int)
            relatedness coefficient and ibs2 (loci with the same genotype)
        """

        # Cache the snp set into the sample and avoid recomputing it every iteration
        if sample1.name not in self.snpSet:
            self.snpSet[sample1.name] = self._createSNPSets(sample1.pileup)

        if sample2.name not in self.snpSet:
            self.snpSet[sample2.name] = self._createSNPSets(sample2.pileup)

        totalCommon = (
            self.snpSet[sample1.name]["total"] & self.snpSet[sample2.name]["total"]
        )

        sharedHets = len(
            self.snpSet[sample1.name]["het"] & self.snpSet[sample2.name]["het"]
        )

        minHet = min(
            len(self.snpSet[sample1.name]["het"] & totalCommon),
            len(self.snpSet[sample2.name]["het"] & totalCommon),
        )

        ibs0 = len(
            self.snpSet[sample1.name]["ref"] & self.snpSet[sample2.name]["alt"]
        ) + len(self.snpSet[sample2.name]["ref"] & self.snpSet[sample1.name]["alt"])

        relatedessCoeff = (sharedHets - 2 * ibs0) / minHet

        # Loci where the two samples have the SAME genotype
        ibs2 = (
            len(self.snpSet[sample1.name]["ref"] & self.snpSet[sample2.name]["ref"])
            + len(self.snpSet[sample2.name]["het"] & self.snpSet[sample1.name]["het"])
            + len(self.snpSet[sample2.name]["alt"] & self.snpSet[sample1.name]["alt"])
        )

        return relatedessCoeff, ibs2

    def _identifyRelation(self, sample1, sample2, relatedness):
        """
        identify the relationship between two samples (twin/replicate/ related/unrelated)
        Args:
            sample1: first sample
            sample2: second sample
            relatedness: relatedness between sample1 and sample2

        """
        if relatedness >= 0.95:
            sample1.replicates += 1
            sample2.replicates += 1
            sample1.replicateSamples.append(sample2.name)
            sample2.replicateSamples.append(sample1.name)
        elif 0.15 <= relatedness < 0.95:
            sample1.familyMembers += 1
            sample2.familyMembers += 1
            sample1.familyMemberSamples.append(sample2.name)
            sample2.familyMemberSamples.append(sample1.name)

    def processRelatedness(self, sample1, sample2):
        """
        compute the relatedness and the relation between two samples
        Args:
            sample1: Sample
                first sample
            sample2: Sample
                second sample

        Returns: Tuple(float, int)
            relatedness coefficient and ibs2 of the pair of samples

        """
        relatedness, ibs2 = self._computeRelatedness(sample1, sample2)
        self._identifyRelation(sample1, sample2, relatedness)
        return relatedness, ibs2
