// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>

namespace OpenMS
{

  //-------------------------------------------------------------
  //Doxygen docu
  //-------------------------------------------------------------

  /**
    @brief Search for cross-linked peptide pairs in tandem MS spectra

    This tool performs a search for cross-links in the given mass spectra.

    It executes the following steps in order:
    <ul>
      <li>Processing of spectra: deisotoping and filtering</li>
      <li>Digesting and preprocessing the protein database, building a peptide pair index dependent on the precursor masses of the MS2 spectra</li>
      <li>Generating theoretical spectra of cross-linked peptides and aligning the experimental spectra against those</li>
      <li>Scoring of cross-link spectrum matches</li>
      <li>Using PeptideIndexer to map the peptides to all possible source proteins</li>
    </ul>

    See below for available parameters and more functionality.

    <h3>Input: MS2 spectra and fasta database of proteins expected to be cross-linked in the sample</h3>
    The spectra should be provided as one PeakMap. If you have multiple files, e.g. for multiple fractions, you should run this tool on each
    file separately.
    The database should be provided as a vector of FASTAEntries containing the target and decoy proteins.

    <h3>Parameters</h3>
    The parameters for fixed and variable modifications refer to additional modifications beside the cross-linker.
    The linker used in the experiment has to be described using the cross-linker specific parameters.
    Only one mass is allowed for a cross-linker that links two peptides, while multiple masses are possible for mono-links of the same cross-linking reagent.
    Mono-links are cross-linkers, that are linked to one peptide by one of their two reactive groups.
    To search for isotopically labeled pairs of cross-linkers see the tool OpenPepXL.
    The parameters -cross_linker:residue1 and -cross_linker:residue2 are used to enumerate the amino acids,
    that each end of the linker can react with. This way any heterobifunctional cross-linker can be defined.
    To define a homobifunctional cross-linker, these two parameters should have the same value.
    The parameter -cross_linker:name is used to solve ambiguities caused by different cross-linkers with the same mass
    after the linking reaction (see section on output for clarification).

    <h3>Output: XL-MS Identifications with scores and linked positions in the proteins</h3>
    The input parameters protein_ids and peptide_ids are filled with XL-MS search parameters and IDs

    <CENTER>
      <table>
          <tr>
              <th ALIGN = "center"> pot. predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> &rarr; OpenPepXLLF &rarr;</td>
              <th ALIGN = "center"> pot. successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
          </tr>
      </table>
    </CENTER>

  */

  class OPENMS_DLLAPI OpenPepXLLFAlgorithm :
   public DefaultParamHandler, public ProgressLogger
  {
public:

    /// Exit codes
    enum ExitCodes
    {
      EXECUTION_OK,
      ILLEGAL_PARAMETERS,
      UNEXPECTED_RESULT,
      INCOMPATIBLE_INPUT_DATA
    };

    /// Default constructor
    OpenPepXLLFAlgorithm();

    /// Default destructor
    ~OpenPepXLLFAlgorithm() override;

    /**
     * @brief Performs the main function of this class, the search for cross-linked peptides

      @param unprocessed_spectra The input PeakMap of experimental spectra
      @param fasta_db The protein database containing targets and decoys
      @param protein_ids A result vector containing search settings. Should contain one PeptideIdentification.
      @param peptide_ids A result vector containing cross-link spectrum matches as PeptideIdentifications and PeptideHits. Should be empty.
      @param all_top_csms A result vector containing cross-link spectrum matches as CrossLinkSpectrumMatches. Should be empty. This is only necessary for writing out xQuest type spectrum files.
      @param spectra A result vector containing the input spectra after preprocessing and filtering. Should be empty. This is only necessary for writing out xQuest type spectrum files.
     */
    ExitCodes run(PeakMap& unprocessed_spectra, std::vector<FASTAFile::FASTAEntry>& fasta_db, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, 
                 std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, PeakMap& spectra);

private:
    void updateMembers_() override;

    String decoy_string_;
    bool decoy_prefix_;

    Int min_precursor_charge_;
    Int max_precursor_charge_;
    double precursor_mass_tolerance_;
    bool precursor_mass_tolerance_unit_ppm_;
    IntList precursor_correction_steps_;

    double fragment_mass_tolerance_;
    double fragment_mass_tolerance_xlinks_;
    bool fragment_mass_tolerance_unit_ppm_;

    StringList cross_link_residue1_;
    StringList cross_link_residue2_;
    double cross_link_mass_;
    DoubleList cross_link_mass_mono_link_;
    String cross_link_name_;

    StringList fixedModNames_;
    StringList varModNames_;
    Size max_variable_mods_per_peptide_;
    Size peptide_min_size_;
    Size missed_cleavages_;
    String enzyme_name_;

    Int number_top_hits_;
    String deisotope_mode_;
    bool use_sequence_tags_;
    Size sequence_tag_min_length_;

    String add_y_ions_;
    String add_b_ions_;
    String add_x_ions_;
    String add_a_ions_;
    String add_c_ions_;
    String add_z_ions_;
    String add_losses_;
  };
}
