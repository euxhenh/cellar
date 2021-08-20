import dash_bootstrap_components as dbc
import dash_html_components as html


lbt_ingest_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Scanpy Ingest Settings"),

        dbc.PopoverBody(
            [
                html.P(
                    [
                        "Details: ",
                        html.A("scanpy.readthedocs.io/en/"
                               "stable/api/scanpy.tl.ingest.html",
                               href="https://scanpy.readthedocs.io/en/"
                               "stable/api/scanpy.tl.ingest.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="lbt-Scanpy-Ingest-settings",
    target="lbt-Scanpy-Ingest-btn",
    trigger="click"
)

lbt_ingest_settings_keys = {

}

lbt_singler_settings = dbc.Popover(
    [
        dbc.PopoverHeader("SingleR Settings"),

        dbc.PopoverBody(
            [
                html.P(
                    [
                        "Details: ",
                        html.A("bioconductor.org/books/release/SingleRBook/",
                               href="http://bioconductor.org/books/release/"
                               "SingleRBook/",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="lbt-SingleR-settings",
    target="lbt-SingleR-btn",
    trigger="click"
)

lbt_singler_settings_keys = {

}

lbt_exact_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Cell ID Based label transfer"),

        dbc.PopoverBody(
            [
                html.P(
                    [
                        "Transfers labels from one dataset to the other "
                        "in case the cell IDs match. Cells without "
                        "a match are assigned a new cluster."
                    ]
                )
            ]
        )
    ],
    id="lbt-exact-settings",
    target="lbt-exact-btn",
    trigger="click"
)

lbt_exact_settings_keys = {

}


lbt_settings = [
    lbt_ingest_settings,
    lbt_singler_settings,
    lbt_exact_settings
]

lbt_settings_keys = {
    'lbt-Scanpy-Ingest': lbt_ingest_settings_keys,
    'lbt-SingleR': lbt_singler_settings_keys,
    'lbt-exact': lbt_exact_settings_keys
}
