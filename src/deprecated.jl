##############################################################################
# deprecated.jl — the names of TensND 0.7 and earlier, renamed in 0.8 to       #
# Julia's convention (a function is lowercase, a proper name included, as      #
# `cholesky` or `schur`). Each still works, with a deprecation warning that    #
# names its successor.                                                         #
##############################################################################

# Kelvin-Mandel matrices and their parameter forms.
@deprecate KM mandel
@deprecate inv_KM from_mandel
@deprecate KM_material mandel_material
@deprecate KM_rotation mandel_rotation
@deprecate KM_from_ti_params mandel_from_ti_params
@deprecate KM_from_ti8_params mandel_from_ti8_params
@deprecate KM_from_ortho_params mandel_from_ortho_params
@deprecate KM_from_cubic_params mandel_from_cubic_params
@deprecate ti_params_from_KM ti_params_from_mandel
@deprecate ti8_params_from_KM ti8_params_from_mandel
@deprecate ortho_params_from_KM ortho_params_from_mandel
@deprecate cubic_params_from_KM cubic_params_from_mandel
@deprecate km_ip_block mandel_ip_block
@deprecate km_op_block mandel_op_block
@deprecate km_ip_embed mandel_ip_embed
@deprecate km_op_embed mandel_op_embed
@deprecate km_from_blocks mandel_from_blocks
@deprecate km_inv_schur mandel_inv_schur
Base.@deprecate_binding KM_IP MANDEL_IP
Base.@deprecate_binding KM_OP MANDEL_OP

# Symmetry predicates.
@deprecate is_ISO isisotropic
@deprecate is_TI istransverselyisotropic
@deprecate is_ORTHO isorthotropic
@deprecate is_CUBIC iscubic

# Transverse isotropy.
@deprecate tens_TI tens_ti
@deprecate tens_TI_eng tens_ti_eng
@deprecate tens_TI_Hoenig tens_ti_hoenig
@deprecate arg_TI ti_args
@deprecate arg_TI_eng ti_args_eng
@deprecate arg_TI_Hoenig ti_args_hoenig
@deprecate fromISO iso_to_ti

# Isotropic projectors and the Walpole basis: the legacy entry points.
@deprecate ISO iso_projectors
function Walpole(n; sym::Bool = false)
    Base.depwarn(
        "`Walpole(n)` is deprecated, use `walpole_basis(n)` " *
            "(or `walpole_basis_sym(n)` for `sym = true`) instead.", :Walpole
    )
    return sym ? walpole_basis_sym(n) : walpole_basis(n)
end
export Walpole

# Special tensors and coordinate systems.
@deprecate LeviCivita levi_civita
@deprecate Christoffel christoffel
@deprecate Lame lame_coefficients
@deprecate getOM position_vector
@deprecate GRAD grad
@deprecate SYMGRAD symgrad
@deprecate DIV divergence
@deprecate LAPLACE laplacian
@deprecate HESS hess
