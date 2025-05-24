module table1

    export do_table1

    using XLSX, DataFrames, PrettyTables, CSV


    function do_table1(m)

        println("main_dir: ", m)

        df = CSV.read(joinpath(m, "data/Private_Education_Participation_Rate_by_School_Level.csv"), DataFrame; header=2)

        #Remove last column
        select!(df, Not("General high school (%)"))

        # List of names to filter out
        names_to_remove = ["Distribution of sample", "Language", "English", "Math", "Social studies, Science", "Essay writing",
                   "Second Language, Chinese, computer, etc.", "Computer", "Second Language, Chinese, etc.",
                   "Music", "Art", "Physical activities", "Hobby and Cultural eudcation", 
                   "Subject: Work-related private tutoring", "Subject: Career Path and Learning Counseling ",
                   "Subjects: General curriculum private education", "Subjects: Arts and Physical Education, hobbies, education private tutoring" ]

        # Loop through the names and filter the rows
        for name in names_to_remove
            filter!(row -> row["Subjects and types"] != name, df)
        end

        # Define the two rows to sum
        name1 = "Textbooks with tutor's visit "
        name2 = "Paid internet and correspondence lecture, etc."

        # Find the rows corresponding to name1 and name2
        row1 = df[df."Subjects and types" .== name1, :]
        row2 = df[df."Subjects and types" .== name2, :]

        # Create an empty row with the same structure as the DataFrame
        empty_row = DataFrame()
        empty_row."Subjects and types" = ["Others"]
        empty_row."Average (%)" = [sum(row1."Average (%)") + sum(row2."Average (%)")]

        # Sum the values for name1 and name2 for each variable and assign to the empty row
        for col in names(df)[2:end]
            empty_row[!, col] = [sum(row1[!, col]) + sum(row2[!, col])]
        end

        # Find the index of the row where "Subjects and types" equals name3
        index_name3 = findfirst(df."Subjects and types" .== "Taking lessons at Private institutes")

        # Append the new row below the row with "Subjects and types" = name3
        if index_name3 !== nothing
            insert!(df, index_name3 + 1, empty_row[1, :])
        else
            append!(df, empty_row)  # This is valid because append! accepts a full DataFrame
        end

        names_to_remove_2 = [name1, name2]

        # Loop through the names and filter the rows
        for name in names_to_remove_2
            filter!(row -> row["Subjects and types"] != name, df)
        end

        # Step 1: Manually format the "Subjects and types" column to match the image
        df[!, :"Subjects and types"] = [
            "Any subject",
            "  A. Main subjects",
            "    a. Individual tutoring",
            "    b. Group tutoring",
            "    c. \\textit{Hagwon}",
            "    d. Others",
            "  B. Art, music, physical activities",
            "    a. Individual tutoring",
            "    b. Group tutoring",
            "    c. \\textit{Hagwon}",
            "    d. Other"
        ]


        # Step 2: Rename the columns for LaTeX table headers
        rename!(df, Dict(
            "Average (%)" => "Average",
            "Elementary school (%)" => "Elementary",
            "Middle school (%)" => "Middle",
            "High school (%)" => "High"
        ))

        # Step 3: Write LaTeX table and footnote to a .tex file
        latex_path = joinpath(m, "output/Table1.tex")
        #"/Users/fra/VS_CCA/Replication/PkgTemplates/output/Table1.tex"

        open(latex_path, "w") do io
        # Write the title as plain LaTeX
        println(io, "\\begin{center}")
        println(io, "\\textbf{Table 1: Private Education Participation Rate (2019, \\%)}")
        println(io, "\\end{center}")
        println(io, "")

        # Write the LaTeX table
        pretty_table(io, df;
            backend = Val(:latex),
            tf = tf_latex_booktabs,
            alignment = :l
        )

        # Write the footnote below the table
        println(io, "\n\\vspace{1em}")
        println(io, "\\begin{minipage}{\\textwidth}")
        println(io, "\\small")
        println(io, "Notes: The main subjects include Korean, English, math, science, second foreign language, writing and computer science. Other activities may include online education programs or home-based sessions with tutors from education companies. Note that the total is not the sum of the lower categories, since a single child may, for example, attend a \\textit{hagwon} and have individual tutoring lessons. Source: \\textit{Private Education Participation Rate by School Level} from \\href{https://kostat.go.kr}{Statistics Korea (2020)}.")
        println(io, "\\end{minipage}")
        end    
    end
end

