extern crate proc_macro;

use proc_macro::*;
use syn::*;
use syn::parse::{Parse, ParseStream};
use quote::ToTokens;

use prime_factorization::Factorization;
use primes::{Sieve, PrimeSet};

struct Number(u128);
struct Range(Number, Number);

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

impl Parse for Range {
    fn parse(input: ParseStream) -> Result<Range> {
        let range = input.parse::<ExprRange>()?;
        let left = range.start
            .ok_or(Error::new(input.span(), "ranges here must be bounded on both sides"))?;
        let left = syn::parse::<Number>(TokenStream::from(left.to_token_stream()))?;
        let right = range.end
            .ok_or(Error::new(input.span(), "ranges here must be bounded on both sides"))?;
        let mut right = syn::parse::<Number>(TokenStream::from(right.to_token_stream()))?;
        if let RangeLimits::Closed(_) = range.limits {
            right.0 += 1;
        }
        Ok(Range(left, right))
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

#[proc_macro]
pub fn primes(tokens: TokenStream) -> TokenStream {
    struct Helper(syn::Ident, Option<syn::Ident>, Range);
    impl Parse for Helper {
        fn parse(input: ParseStream) -> Result<Helper> {
            let callback = input.parse::<syn::Ident>()?;
            let marker = if input.peek2(syn::Ident) {
                input.parse::<Token![,]>()?;
                Some(input.parse::<syn::Ident>()?)
            } else {
                None
            };
            input.parse::<Token![,]>()?;
            let range = input.parse::<Range>()?;
            Ok(Helper(callback, marker, range))
        }
    }
    let Helper(callback, marker, Range(start, end)) = parse_macro_input!(tokens as Helper);
    let mut args = Vec::<TokenTree>::new();
    if let Some(marker) = marker {
        TokenStream::from(marker.to_token_stream())
            .into_iter()
            .for_each(|x| args.push(x));
        args.push(TokenTree::Punct(Punct::new(',', Spacing::Alone)));
    }
    Sieve::new()
        .iter()
        .skip_while(|x| x < &(start.0 as u64))
        .take_while(|x| x < &(end.0 as u64))
        .for_each(|x| {
            args.push(TokenTree::Literal(Literal::u128_unsuffixed(x as u128)));
            args.push(TokenTree::Punct(Punct::new(',', Spacing::Alone)));
        });
    TokenStream::from_iter(vec![
        TokenTree::Ident(proc_macro::Ident::new(&callback.to_string(), Span::call_site())),
        TokenTree::Punct(Punct::new('!', Spacing::Alone)),
        TokenTree::Group(Group::new(
            Delimiter::Parenthesis,
            TokenStream::from_iter(args),
        )),
        TokenTree::Punct(Punct::new(';', Spacing::Alone))
    ])
}

