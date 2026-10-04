$ErrorActionPreference = 'Stop'

$input = Join-Path $PSScriptRoot 'lactose_synthase_glucose_control_1NF5_original.pdb'
$receptorOutput = Join-Path $PSScriptRoot 'lactose_synthase_1NF5_AB_receptor_no_glucose.pdb'
$referenceOutput = Join-Path $PSScriptRoot 'lactose_synthase_1NF5_bound_beta_D_glucose_reference.pdb'

$lines = Get-Content -LiteralPath $input
$receptor = New-Object System.Collections.Generic.List[string]
$glucose = New-Object System.Collections.Generic.List[string]

foreach ($line in $lines) {
    if ($line.Length -lt 27) { continue }
    $record = $line.Substring(0, 6).Trim()
    $residue = $line.Substring(17, 3).Trim()
    $chain = $line.Substring(21, 1)
    $residueNumber = $line.Substring(22, 4).Trim()

    # One experimental lactose-synthase complex: alpha-lactalbumin A + beta4Gal-T1 B.
    if ($record -eq 'ATOM' -and $chain -in @('A', 'B')) {
        $receptor.Add($line)
    }

    # Preserve the alpha-lactalbumin structural calcium, CA A 124.
    if ($record -eq 'HETATM' -and $residue -eq 'CA' -and $chain -eq 'A' -and $residueNumber -eq '124') {
        $receptor.Add($line)
    }

    # Experimental bound glucose: beta-D-glucopyranose, BGC B 403.
    if ($record -eq 'HETATM' -and $residue -eq 'BGC' -and $chain -eq 'B' -and $residueNumber -eq '403') {
        $glucose.Add($line)
    }
}

$receptor.Add('TER')
$receptor.Add('END')
$glucose.Add('END')

Set-Content -LiteralPath $receptorOutput -Value $receptor -Encoding ascii
Set-Content -LiteralPath $referenceOutput -Value $glucose -Encoding ascii

Write-Output "Wrote receptor: $receptorOutput"
Write-Output "Wrote experimental glucose reference: $referenceOutput"
Write-Output "Removed from receptor: BGC chain B residue 403 (12 atoms)"
