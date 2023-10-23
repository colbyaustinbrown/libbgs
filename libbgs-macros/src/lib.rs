extern crate proc_macro;

use proc_macro::*;
use syn::*;
use syn::parse::{Parse, ParseStream};

use prime_factorization::Factorization;

struct Number(u128);

impl Parse for Number {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(token::Brace) {
            let expr;
            braced!(expr in input);
            let bin = expr.parse::<ExprBinary>()?;
            let mut left = *bin.left;
            if let Expr::Group(ExprGroup {
                expr,
                ..
            }) = left {
                left = *expr;
            };
            let Expr::Lit(ExprLit {
                lit: Lit::Int(a),
                attrs: _,
            }) = left else {
                return Err(Error::new(expr.span(), format!("expected an unsigned integer literal, instead saw {:?}", left)));
            };
            let mut right = *bin.right;
            if let Expr::Group(ExprGroup {
                expr,
                ..
            }) = right {
                right = *expr;
            };
            let Expr::Lit(ExprLit {
                lit: Lit::Int(b),
                attrs: _,
            }) = right else {
                return Err(Error::new(expr.span(), "expected an unsigned integer literal b"));
            };
            let a = LitInt::from(a).base10_parse::<u128>()?;
            let b = LitInt::from(b).base10_parse::<u128>()?;
            match bin.op {
                BinOp::Add(_) => Ok(Number(a + b)),
                BinOp::Sub(_) => Ok(Number(a - b)),
                _ => {
                    return Err(Error::new(expr.span(), "only addition or subtraction allowed here"));
                }
            }
        } else {
            let num = input.parse::<LitInt>()?;
            Ok(Number(num.base10_parse::<u128>()?))
        }
    }
}

#[proc_macro]
pub fn make_factor(tokens: TokenStream) -> TokenStream {
    let mut res = Vec::<TokenTree>::new();
    let Number(n) = syn::parse::<Number>(tokens).unwrap();

    res.push(TokenTree::Punct(Punct::new('&', Spacing::Alone)));
    let mut entries = Vec::<TokenTree>::new();
    for (p, t) in Factorization::run(n).prime_factor_repr() {
        entries.push(TokenTree::Group(Group::new(
            Delimiter::Parenthesis,
            TokenStream::from_iter(vec![
                TokenTree::Literal(Literal::u128_unsuffixed(p)),
                TokenTree::Punct(Punct::new(',', Spacing::Alone)),
                TokenTree::Literal(Literal::u128_unsuffixed(t as u128)),
            ]),
        )));
        entries.push(TokenTree::Punct(Punct::new(',', Spacing::Alone)));
    }
    res.push(TokenTree::Group(Group::new(
        Delimiter::Bracket,
        TokenStream::from_iter(entries),
    )));

    TokenStream::from_iter(res)
}
